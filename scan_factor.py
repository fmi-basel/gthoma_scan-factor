# Copyright 2024 Friedrich Miescher Institute for Biomedical Research
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Authors: Richard Bunker, Georg Kempf, Friedrich Miescher Institute for Biomedical Research

import os
import time
import traceback

import pandas as pd
from matplotlib.figure import Figure
from pymol import cmd

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
_TRUE_STRINGS = {"1", "true", "yes", "y", "on"}
_ATOM_NAME_ALIASES = {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3"}
_MIN_ATOM_PAIRS = 3

RESULT_COLUMNS = [
    "Atom clashes Protein/Protein",
    "Residue clashes Protein/Protein",
    "Atom clashes Protein/DNA",
    "Residue clashes Protein/DNA",
    "Residue clashes overall",
    "Atom clashes overall",
    "RMSD",
]

# Hidden helper selections (all deleted at the end via "_sf_*")
_SEL_PROBE_PROT = "_sf_probe_prot"
_SEL_TARGET_PROT = "_sf_target_prot"
_SEL_TARGET_DNA = "_sf_target_dna"


# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def _to_bool(value):
    """PyMOL passes command-line arguments as strings ('False' would be truthy)."""
    if isinstance(value, str):
        return value.strip().lower() in _TRUE_STRINGS
    return bool(value)


def _to_optional_int(value):
    if value is None or (isinstance(value, str) and value.strip().lower() in ("", "none")):
        return None
    return int(value)


def _fmt_time(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:d}:{m:02d}:{s:02d}" if h else f"{m:d}:{s:02d}"


# --------------------------------------------------------------------------- #
# Progress dialog
# --------------------------------------------------------------------------- #
class _Progress:
    """
    Progress dialog for long PyMOL loops.
    Uses a Qt QProgressDialog when the PyMOL GUI is available and we are on the
    GUI thread; otherwise falls back to console output (e.g. `pymol -c`).
    """

    def __init__(self, total, title="scan_factor", enabled=True):
        self.total = max(int(total), 1)
        self.t0 = time.time()
        self.dialog = None
        self._app = None
        self._last_pct = -1
        if enabled:
            self._init_qt(title)

    def _init_qt(self, title):
        try:
            from pymol.Qt import QtWidgets, QtCore
        except ImportError:
            return
        app = QtWidgets.QApplication.instance()
        # Qt widgets must only be created on the GUI thread
        if app is None or QtCore.QThread.currentThread() != app.thread():
            return
        try:
            dlg = QtWidgets.QProgressDialog("Preparing…", "Cancel", 0, self.total,
                                            app.activeWindow())
            dlg.setWindowTitle(title)
            # Qt5 / Qt6 compatible enum lookup
            modal = getattr(QtCore.Qt, "WindowModal", None) or QtCore.Qt.WindowModality.WindowModal
            dlg.setWindowModality(modal)
            dlg.setMinimumDuration(0)
            dlg.setAutoClose(False)
            dlg.setAutoReset(False)
            dlg.setMinimumWidth(420)
            dlg.setValue(0)
            dlg.show()
            app.processEvents()
            self.dialog, self._app = dlg, app
        except Exception as exc:  # never let the GUI break the scan
            print(f"Progress dialog unavailable ({exc}); using console output.")
            self.dialog = None

    def update(self, done, text=""):
        elapsed = time.time() - self.t0
        remaining = elapsed / done * (self.total - done) if done else 0.0
        status = (f"{done}/{self.total} positions  |  elapsed {_fmt_time(elapsed)}"
                  f"  |  remaining ~{_fmt_time(remaining)}")
        if self.dialog is not None:
            self.dialog.setLabelText(f"{text}\n{status}" if text else status)
            self.dialog.setValue(done)
            self._app.processEvents()  # keeps dialog and 3D viewer responsive
        else:
            pct = int(100 * done / self.total)
            if pct // 5 != self._last_pct // 5 or done == self.total:  # every 5 %
                self._last_pct = pct
                print(f"[{pct:3d}%] {status}")

    @property
    def cancelled(self):
        return self.dialog is not None and self.dialog.wasCanceled()

    def close(self):
        if self.dialog is not None:
            self.dialog.close()
            self.dialog.deleteLater()
            self._app.processEvents()
            self.dialog = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# --------------------------------------------------------------------------- #
# PyMOL helpers
# --------------------------------------------------------------------------- #
def _residue_numbers(obj, chain):
    """Ascending unique residue numbers of nucleic-acid residues in obj/chain."""
    out = set()
    cmd.iterate(f"({obj}) and polymer.nucleic and chain {chain}",
                "out.add(resv)", space={"out": out})
    return sorted(out)


def _backbone_atoms(obj):
    """Return {(chain, resv): {atom_name: index}} for nucleic-acid backbone atoms."""
    rows = []
    cmd.iterate(f"({obj}) and polymer.nucleic and backbone",
                "rows.append((chain, resv, name, alt, index))",
                space={"rows": rows})
    atoms = {}
    for chain, resv, name, alt, index in rows:
        if alt not in ("", "A"):  # keep only one alternate conformation
            continue
        name = _ATOM_NAME_ALIASES.get(name, name)
        atoms.setdefault((chain, resv), {}).setdefault(name, index)
    return atoms


def _count_residues(selection):
    res = set()
    cmd.iterate(selection, "res.add((segi, chain, resi))", space={"res": res})
    return len(res)


def _build_atom_pairs(residue_pairs, probe, target, probe_bb, target_bb, label):
    """
    Build an explicit, ordered list of atom selections for cmd.pair_fit
    (probe_atom, target_atom, probe_atom, target_atom, ...), using only
    backbone atoms present in both residues of each pair.
    """
    selections, probe_indices = [], []
    for p_key, t_key in residue_pairs:
        p_atoms = probe_bb.get(p_key, {})
        t_atoms = target_bb.get(t_key, {})
        unmatched = set(p_atoms) ^ set(t_atoms)
        if unmatched:
            print(f"Warning ({label}): omitting unmatched backbone atoms {sorted(unmatched)} "
                  f"(probe {p_key[0]}/{p_key[1]} vs target {t_key[0]}/{t_key[1]})")
        for name in sorted(set(p_atoms) & set(t_atoms)):
            selections.append(f"({probe}) and index {p_atoms[name]}")
            selections.append(f"({target}) and index {t_atoms[name]}")
            probe_indices.append(p_atoms[name])
    return selections, probe_indices


# --------------------------------------------------------------------------- #
# Scan
# --------------------------------------------------------------------------- #
def _scan(target, probe, first_base_probe, probe_center, probe_length, cutoff,
          clash_keep, write_models, align_both_strands, show_progress=True):

    # ---- argument checks ---------------------------------------------------
    objects = cmd.get_names("objects")
    for obj in (target, probe):
        if obj not in objects:
            raise ValueError(f'Object "{obj}" not found.')
    if probe_center and first_base_probe is not None:
        raise ValueError("Arguments first_base_probe and probe_center cannot be used at the same time.")
    if probe_length < 1:
        raise ValueError("probe_length must be >= 1.")

    # ---- clean up (only the two objects involved) -------------------------
    print("Removing solvent and hydrogens")
    cmd.remove(f"(({target}) or ({probe})) and (solvent or hydro)")

    # ---- DNA chains --------------------------------------------------------
    target_chains = cmd.get_chains(f"({target}) and polymer.nucleic")
    probe_chains = cmd.get_chains(f"({probe}) and polymer.nucleic")
    if len(target_chains) != 2:
        raise ValueError(f'Target "{target}" must contain exactly two DNA strands '
                         f'(found {len(target_chains)}).')
    if len(probe_chains) != 2:
        raise ValueError(f'Probe "{probe}" must contain exactly two DNA strands '
                         f'(found {len(probe_chains)}).')
    t1, t2 = target_chains
    p1, p2 = probe_chains

    # ---- probe window ------------------------------------------------------
    # Strand 1 5'->3'; strand 2 reversed so that index i on both lists is a base pair
    p1_res = _residue_numbers(probe, p1)
    p2_res = _residue_numbers(probe, p2)[::-1]
    if len(p1_res) != len(p2_res):
        raise ValueError("The two probe DNA strands have different lengths.")
    if probe_length > len(p1_res):
        raise ValueError(f"probe_length ({probe_length}) exceeds probe DNA length ({len(p1_res)}).")

    if probe_center:
        start = (len(p1_res) - probe_length) // 2
    elif first_base_probe is not None:
        if first_base_probe not in p1_res:
            raise ValueError(f"first_base_probe {first_base_probe} not found in probe chain {p1}.")
        start = p1_res.index(first_base_probe)
    else:
        start = 0
    if start + probe_length > len(p1_res):
        raise ValueError("Probe window extends beyond the end of the probe DNA.")

    probe_win_1 = p1_res[start:start + probe_length]
    probe_win_2 = p2_res[start:start + probe_length]

    # ---- target residues ---------------------------------------------------
    t1_res = _residue_numbers(target, t1)
    t2_res = _residue_numbers(target, t2)
    if align_both_strands and len(t1_res) != len(t2_res):
        raise ValueError("Target DNA strands have different lengths; cannot align both strands.")

    print(f"Probe residues used in alignment:\n"
          f"  Chain {p1}: {probe_win_1[0]}-{probe_win_1[-1]}\n"
          f"  Chain {p2}: {min(probe_win_2)}-{max(probe_win_2)}\n"
          f"Probe length: {probe_length}\n"
          f"Target length: {len(t1_res)} ({t1}) / {len(t2_res)} ({t2})")

    # ---- build scan positions ---------------------------------------------
    jobs = []
    if align_both_strands:
        t2_rev = t2_res[::-1]
        for i in range(len(t1_res) - probe_length + 1):
            w1 = t1_res[i:i + probe_length]
            w2 = t2_rev[i:i + probe_length]
            pairs = [((p1, a), (t1, b)) for a, b in zip(probe_win_1, w1)]
            pairs += [((p2, a), (t2, b)) for a, b in zip(probe_win_2, w2)]
            jobs.append({"chain": f"{t1}_{t2}", "start": w1[0], "end": w1[-1],
                         "pose": f"{t1}_{w1[0]}_{w1[-1]}_{t2}_{w2[-1]}_{w2[0]}",
                         "pairs": pairs})
    else:
        # Probe strand 1 is scanned along each target strand (both orientations)
        for t_chain, t_res in ((t1, t1_res), (t2, t2_res)):
            for i in range(len(t_res) - probe_length + 1):
                w = t_res[i:i + probe_length]
                pairs = [((p1, a), (t_chain, b)) for a, b in zip(probe_win_1, w)]
                jobs.append({"chain": t_chain, "start": w[0], "end": w[-1],
                             "pose": f"{t_chain}_{w[0]}_{w[-1]}", "pairs": pairs})

    # ---- pre-compute atoms and helper selections ---------------------------
    probe_bb = _backbone_atoms(probe)
    target_bb = _backbone_atoms(target)
    cmd.select(_SEL_PROBE_PROT, f"({probe}) and polymer.protein", enable=0)
    cmd.select(_SEL_TARGET_PROT, f"({target}) and polymer.protein", enable=0)
    cmd.select(_SEL_TARGET_DNA, f"({target}) and polymer.nucleic", enable=0)

    group_name = f"clash_keep_{clash_keep}_probe_length_{probe_length}"
    model_dir = f"{target}_{probe}_models"
    if write_models:
        os.makedirs(model_dir, exist_ok=True)

    sel_pp = f"{_SEL_PROBE_PROT} within {cutoff} of {_SEL_TARGET_PROT}"
    sel_pd = f"{_SEL_PROBE_PROT} within {cutoff} of {_SEL_TARGET_DNA}"

    # ---- main loop ---------------------------------------------------------
    results = []
    cancelled = False
    print(f"Scanning {len(jobs)} positions")

    with _Progress(len(jobs), f"scan_factor: {probe} → {target}", show_progress) as progress:
        try:
            for n, job in enumerate(jobs):
                label = f"DNA chain {job['chain']} bases {job['start']}-{job['end']}"

                if progress.cancelled:
                    cancelled = True
                    print(f"Scan cancelled by user after {n}/{len(jobs)} positions.")
                    break
                progress.update(n, f"Aligning {probe} on {label}")

                row = {"chain": job["chain"], "start": job["start"], "end": job["end"]}

                atom_sels, probe_indices = _build_atom_pairs(
                    job["pairs"], probe, target, probe_bb, target_bb, label)

                rmsd = None
                if len(atom_sels) // 2 >= _MIN_ATOM_PAIRS:
                    try:
                        rmsd = cmd.pair_fit(*atom_sels, quiet=1)
                    except Exception as exc:
                        print(f"{label}: pair_fit error: {exc}")

                if rmsd is None or rmsd < 0:
                    print(f"{label}: {probe} superposition failed, clashes not calculated – check for "
                          f"missing atoms, alternate conformations or unusual bases at these positions")
                    row.update({col: None for col in RESULT_COLUMNS})
                    results.append(row)
                    continue

                # Clashes of probe protein with target protein / target DNA
                pp_atoms = cmd.count_atoms(sel_pp)
                pd_atoms = cmd.count_atoms(sel_pd)
                pp_res = _count_residues(sel_pp)
                pd_res = _count_residues(sel_pd)

                # Write full models (target + superimposed probe), e.g. for Rosetta scoring
                if write_models:
                    cmd.save(os.path.join(model_dir, f"{target}_{probe}_{job['pose']}.pdb"),
                             f"({target}) or ({probe})")

                # Keep poses with protein/protein residue clashes <= clash_keep
                if pp_res <= clash_keep:
                    obj_name = cmd.get_legal_name(f"probe_{job['pose']}")
                    idx = "+".join(map(str, probe_indices))
                    cmd.create(obj_name,
                               f"({probe}) and (polymer.protein or (byres (({probe}) and index {idx})))")
                    cmd.group(group_name, obj_name)

                print(f"{label}: RMSD {rmsd:.2f} Å | protein/protein {pp_atoms} atoms in {pp_res} residues"
                      f" | protein/DNA {pd_atoms} atoms in {pd_res} residues")

                row.update({
                    "Atom clashes Protein/Protein": pp_atoms,
                    "Residue clashes Protein/Protein": pp_res,
                    "Atom clashes Protein/DNA": pd_atoms,
                    "Residue clashes Protein/DNA": pd_res,
                    "Residue clashes overall": pp_res + pd_res,
                    "Atom clashes overall": pp_atoms + pd_atoms,
                    "RMSD": rmsd,
                })
                results.append(row)
            else:
                progress.update(len(jobs), "Done")
        finally:
            cmd.delete("_sf_*")

    df = pd.DataFrame(results, columns=["chain", "start", "end"] + RESULT_COLUMNS)
    return df, cancelled


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def _plot(df, target, probe, probe_length):
    df = df.copy()
    for col in ["start"] + RESULT_COLUMNS:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["start", "RMSD"])
    if df.empty:
        print("No successful superpositions – nothing to plot.")
        return

    for chain, group in df.groupby("chain"):
        group = group.sort_values("start")
        fig = Figure(figsize=(10, 5))
        ax = fig.subplots()
        ax.plot(group["start"], group["Residue clashes Protein/Protein"], label="Protein/Protein")
        ax.plot(group["start"], group["Residue clashes Protein/DNA"], label="Protein/DNA")
        ax.plot(group["start"], group["Residue clashes overall"], label="Overall")
        ax.set_xlabel("Start base (target)")
        ax.set_ylabel("Residue clashes")
        ax.set_title(f"Residue clashes for chain {chain}")
        ax.legend(loc="upper left")

        ax2 = ax.twinx()
        ax2.plot(group["start"], group["RMSD"], color="red", linestyle="--", label="RMSD")
        ax2.set_ylabel("RMSD (Å)")
        ax2.legend(loc="upper right")

        fig.tight_layout()
        out = f"{target}_{probe}_residue_clashes_combined_{chain}_probe_length_{probe_length}.png"
        fig.savefig(out, dpi=150)
        print(f"Saved plot: {out}")


# --------------------------------------------------------------------------- #
# PyMOL command
# --------------------------------------------------------------------------- #
@cmd.extend
def scan_factor(target, probe, first_base_probe=None, probe_center=False, probe_length=2,
                cutoff=2, clash_keep=100000, write_models=False, align_both_strands=False,
                only_plot=False, show_progress=True):
    """
    DESCRIPTION

        Iteratively align the DNA of a probe object (e.g. TF bound to DNA motif) onto
        the DNA of a target object (e.g. nucleosome) and calculate the clashes between
        protein atoms of the probe and DNA and/or protein atoms of the target.

    USAGE

        scan_factor target, probe [, first_base_probe [, probe_center [, probe_length
                    [, cutoff [, clash_keep [, write_models [, align_both_strands
                    [, only_plot [, show_progress ]]]]]]]]]

    ARGUMENTS

        target = object
        probe = object
        first_base_probe = int: residue number in probe DNA (strand 1) where the
                           superposition selection starts.
        probe_center = bool: centre the superposition window on the probe DNA
                       (default False). Cannot be combined with first_base_probe.
        probe_length = int: number of bases used for superposition (default 2)
        cutoff = float: distance (Å) below which atoms are considered clashing (default 2)
        clash_keep = int: keep poses with <= this many protein/protein residue clashes
                     (default 100000)
        write_models = bool: write target + superimposed probe models to
                       ./<target>_<probe>_models (default False)
        align_both_strands = bool: superimpose using both DNA strands (default False)
        only_plot = bool: only generate plots from an existing CSV (default False)
        show_progress = bool: show a progress dialog with a Cancel button
                        (falls back to console output without GUI; default True)

    EXAMPLE

        fetch 6t90, async=0
        create NCP, 6t90 and c. A+B+C+D+E+F+G+H+I+J
        create OCTSOX, 6t90 and (c. K+L or (c. I and i. 3-24) or (c. J and i. 124-145))
        scan_factor NCP, OCTSOX, probe_center=True
    """
    try:
        probe_center = _to_bool(probe_center)
        write_models = _to_bool(write_models)
        align_both_strands = _to_bool(align_both_strands)
        only_plot = _to_bool(only_plot)
        show_progress = _to_bool(show_progress)
        first_base_probe = _to_optional_int(first_base_probe)
        probe_length = int(probe_length)
        cutoff = float(cutoff)
        clash_keep = int(clash_keep)

        csv_out = f"{target}_{probe}_probe_length_{probe_length}_clashes.csv"

        if only_plot:
            if not os.path.exists(csv_out):
                raise FileNotFoundError(f"CSV not found: {csv_out}")
            df = pd.read_csv(csv_out)
        else:
            df, cancelled = _scan(target, probe, first_base_probe, probe_center, probe_length,
                                  cutoff, clash_keep, write_models, align_both_strands,
                                  show_progress)
            if cancelled:
                # don't overwrite results from a previous complete run
                csv_out = csv_out.replace("_clashes.csv", "_clashes_partial.csv")
            print(df.to_string())
            df.to_csv(csv_out, index=False)
            print(f"Results written to {csv_out}")

        _plot(df, target, probe, probe_length)

    except Exception as exc:
        print(f"scan_factor failed: {exc}")
        traceback.print_exc()