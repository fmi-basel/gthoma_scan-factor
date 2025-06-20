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


import sys, os
import traceback
from pymol import cmd, stored, editing
import pandas as pd
import matplotlib.pyplot as plt

@cmd.extend
def scan_factor(target, probe, first_base_probe=None, probe_center=False, probe_length=2, cutoff=2, clash_keep=100000, write_models=False, align_both_strands=False, only_plot=False):
  """

  DESCRIPTION

    Iteratively align the DNA of a probe object (e.g. TF bound to DNA motif) onto the DNA of a target object (e.g. nucleosome) and calculate the clashes between protein atoms of probe and DNA and/or protein atoms of target.

  USAGE

      scan_factor target, probe [, first_base_probe [, probe_center [, probe_length [, cutoff [, clash_keep [, write_models [, align_both_strands [, only_plot ]]]]]]]]

  ARGUMENTS

      target = object
      probe = object
      first_base_probe = int: residue number in probe DNA where to start selection for superposition.
      probe_center = bool: Use half-length of probe DNA - half of probe_length (rounded) as starting point for superposition selection (default False). Cannot be used together with first_base_probe.
      probe_length = int: number of bases on probe to use for superposition (default 2)
      cutoff = int: distance in Ångstrom below which atoms are considered clashing (default 2)
      clash_keep = int: only keep models with fewer than the specified number of residue clashes (default 100000)
      write_models = bool: Write models below clash_keep threshold to current folder
      align_both_strands=bool: Do the superposition with both DNA strands
      only_plot = bool: Only generate plots based on existing results (csv) in the current folder
      
  EXAMPLE

      fetch 6t90, async=0
      create NCP, 6t90 and c. A+B+C+D+E+F+G+H+I+J
      create OCTSOX, 6t90 and (c. K+L or (c. I and i. 3-24) or (c. J and i. 124-145))
      scan_factor NCP, OCTSOX, probe_center=True

  """
  try:

    def get_key_by_value(dict_, search_value):
      """Get dict key by value"""
      matching_keys = []
      for key, value in dict_.items():
        if value == search_value:
            matching_keys.append(key)
      if len(matching_keys) > 1:
        print(matching_keys)
        raise ValueError("More than one matching values, expected values to be unique.")
      else:
        return matching_keys[0]

    results = []

    #Remove solvent, like water
    cmd.remove("solvent")
    if write_models:
      dir_name = (f"{target}_{probe}_models")
      if not os.path.exists(dir_name):
        os.mkdir("%s" % (dir_name) )

    #Check if DNA is present in probe
    if not cmd.select(f"polymer.nucleic and {target}"):
      print(f"The target object \"{target}\" must contain contain DNA – are you sure this is a target model?\n")

    #Check if DNA is present in target
    if not cmd.select(f"polymer.nucleic and \"{probe}\""):
      print(f"The probe object \"{probe}\" must contain at least two DNA bases for superposition,"
            f" the first two of which should correspond "
            f"to the position of the central two bases of the factor's recognition sequence when bound to DNA\n")

    # Determine probe DNA chain ID
    cmd.select("probeDNA", f"polymer.nucleic and {probe}")
    probe_chains = cmd.get_chains(f"polymer.nucleic and {probe}")
    probe_chain_1 = probe_chains[0]
    probe_chain_2 = probe_chains[1]
    probe_length = int(probe_length)
    #Get chains
    probe_object_chain_1 = cmd.get_model(f"{probe} and chain {probe_chain_1}")
    probe_object_chain_2 = cmd.get_model(f"{probe} and chain {probe_chain_2}")
    #Get residues and store in list
    probe_residues_chain_1 = sorted(list(set([int(at.resi) for at in probe_object_chain_1.atom])))
    probe_residues_chain_2 = sorted(list(set([int(at.resi) for at in probe_object_chain_2.atom])), reverse=True)
    #Map residue numbers to an index starting at 0
    probe_residues_chain_1_mapping = {i: int(resi) for (i, resi) in enumerate(probe_residues_chain_1)}
    probe_residues_chain_2_mapping = {i: int(resi) for (i, resi) in enumerate(probe_residues_chain_2)}
    #Debugging
    #print(probe_residues_chain_1_mapping)
    #print(probe_residues_chain_2_mapping)
    #print(probe_residues_chain_1_mapping)
    #print(probe_residues_chain_2_mapping)
    assert(len(probe_residues_chain_1_mapping) == len(probe_residues_chain_2_mapping))

    #Get starting and end bases of probe
    if first_base_probe and probe_center:
        raise ValueError("Arguments first_base_probe and probe center cannot be used at the same time")
    if probe_center:
      len_center = int(len(probe_residues_chain_1) / 2 ) - int(probe_length / 2)
      first_base_probe_chain_1 = probe_residues_chain_1_mapping[len_center - 1]
      key_first_base = get_key_by_value(probe_residues_chain_1_mapping, first_base_probe_chain_1)
      last_base_probe_chain_2 = probe_residues_chain_2_mapping[key_first_base]
    elif first_base_probe:
      first_base_probe_chain_1 = int(first_base_probe)
      key_first_base = get_key_by_value(probe_residues_chain_1_mapping, first_base_probe_chain_1)
      last_base_probe_chain_2 = probe_residues_chain_2_mapping[key_first_base]
    else:
      key_first_base = 0
      first_base_probe_chain_1 = probe_residues_chain_1_mapping[key_first_base]
      last_base_probe_chain_2 = probe_residues_chain_2_mapping[key_first_base]
    last_base_probe_chain_1_index = key_first_base + probe_length - 1
    last_base_probe_chain_1 = probe_residues_chain_1_mapping[last_base_probe_chain_1_index]
    first_base_probe_chain_2 = probe_residues_chain_2_mapping[last_base_probe_chain_1_index]

    # Remove hydrogens
    print ("Removing hydrogens\n")
    cmd.remove("hydrogens")
    
    # Determine target DNA chain IDs
    target_chains = cmd.get_chains(f"polymer.nucleic and {target}")
    if len(target_chains) < 2:
      raise ValueError("Expected two DNA strands in target")
    elif len(target_chains) > 2:
      raise ValueError("Expected not more than two DNA strands in target")


    target_chain_1 = target_chains[0]
    target_chain_2 = target_chains[1]

    target_atoms_chain_1 = cmd.get_model(f"{target} and chain {target_chain_1}")
    target_atoms_chain_2 = cmd.get_model(f"{target} and chain {target_chain_2}")

    target_residues_chain_1 = sorted(list(set([int(at.resi) for at in target_atoms_chain_1.atom])))
    target_residues_chain_2 = sorted(list(set([int(at.resi) for at in target_atoms_chain_2.atom])), reverse=True)
    target_residues_chain_1_mapping = {i: int(resi) for (i, resi) in enumerate(target_residues_chain_1)}
    target_residues_chain_2_mapping = {i: int(resi) for (i, resi) in enumerate(target_residues_chain_2)}
    assert(len(target_residues_chain_1_mapping) == len(target_residues_chain_2_mapping))
    #Debugging
      #print(f"Target residues chain 1: {target_chain_1}")
    #print(target_residues_chain_1)
    #print("Target residue mapping chain 1:")
    #print(target_residues_chain_1_mapping)
    #print("Target residue mapping chain 2:")
    #print(target_residues_chain_2_mapping)

    first_base_target_chain_1 = target_residues_chain_1_mapping[0]
    first_base_target_chain_2 = target_residues_chain_2_mapping[0]
    target_length = len(target_residues_chain_1_mapping)

    print(f"Probe residues used in alignment:\nChain1:{first_base_probe_chain_1}-{last_base_probe_chain_1}\nChain2:{first_base_probe_chain_2}-{last_base_probe_chain_2}\nProbe length: {probe_length}\nTarget length: {target_length}")

    if align_both_strands:
       target_chains = [target_chain_1]
    else:
       target_chains = [target_chain_1, target_chain_2]

    if not only_plot:
      # Iterate over the residues
      print("Iterating over target DNA strand")
      for chain_idx, target_chain in enumerate(target_chains):
        for first_base_target_chain_1_index in range(0, target_length - probe_length - 1):
          print(f"Current base: {first_base_target_chain_1_index}")

          #Forward DNA
          last_base_target_chain_1_index = (first_base_target_chain_1_index + int(probe_length) - 1)
          first_base_target_chain_1 = target_residues_chain_1_mapping[first_base_target_chain_1_index]
          last_base_target_chain_1 = target_residues_chain_1_mapping[last_base_target_chain_1_index]

          #Reverse DNA
          last_base_target_chain_2 = target_residues_chain_2_mapping[first_base_target_chain_1_index]
          first_base_target_chain_2 = target_residues_chain_2_mapping[last_base_target_chain_1_index]

          # Add \ to negative residue numbers
          if first_base_target_chain_1 < 0:
             first_base_target_chain_1 = f'\{first_base_target_chain_1}'
          if last_base_target_chain_1 < 0:
             last_base_target_chain_1 = f'\{last_base_target_chain_1}'
          if first_base_target_chain_2 < 0:
             first_base_target_chain_2 = f'\{first_base_target_chain_2}'
          if last_base_target_chain_2 < 0:
             last_base_target_chain_2 = f'\{last_base_target_chain_2}'
          if first_base_probe_chain_1 < 0:
              first_base_probe_chain_1 = f'\{first_base_probe_chain_1}'
          if last_base_probe_chain_1 < 0:
             last_base_probe_chain_1 = f'\{last_base_probe_chain_1}'   
          if first_base_probe_chain_2 < 0:
              first_base_probe_chain_2 = f'\{first_base_probe_chain_2}'
          if last_base_probe_chain_2 < 0:
             last_base_probe_chain_2 = f'\{last_base_probe_chain_2}' 
          
          #Debugging
          #print(f"first_base of target forward strand: {first_base_target_chain_1}")
          #print(f"last_base of target forward strand: {last_base_target_chain_1}")
          #print(f"first base of target reverse strand: {first_base_target_chain_2}")
          #print(f"last base of target reverse strand: {last_base_target_chain_2}")
          #print(first_base_target_chain_1)
          #print(last_base_target_chain_1)
          
          #Select residue ranges for current chain
          if not align_both_strands:
            print(chain_idx)
            if chain_idx == 0:
               target_chain = target_chain_1
               probe_chain = probe_chain_1
               first_base_target_chain = first_base_target_chain_1
               last_base_target_chain = last_base_target_chain_1
               first_base_probe_chain = first_base_probe_chain_1
               last_base_probe_chain = last_base_probe_chain_1
            elif chain_idx == 1:
              target_chain = target_chain_2
              probe_chain = probe_chain_1
              first_base_target_chain = first_base_target_chain_2
              last_base_target_chain = last_base_target_chain_2
              first_base_probe_chain = first_base_probe_chain_1
              last_base_probe_chain = last_base_probe_chain_1
            else:
               raise ValueError('Only two chains allowed')

          print(probe_chain_1)
          if align_both_strands:
            moving_selection = f"{probe} and ((chain {probe_chain_1} and resi {first_base_probe_chain_1}-{last_base_probe_chain_1}) or (chain {probe_chain_2} and resi {first_base_probe_chain_2}-{last_base_probe_chain_2})) and backbone"
          else:
            moving_selection = f"{probe} and chain {probe_chain} and resi {first_base_probe_chain}-{last_base_probe_chain} and backbone"
          #Debugging
          #print("Moving selection")
          #print(moving_selection)
          if align_both_strands:
            target_selection = f"{target} and ((chain {target_chain_1} and resi {first_base_target_chain_1}-{last_base_target_chain_1}) or (chain {target_chain_2} and resi {first_base_target_chain_2}-{last_base_target_chain_2})) and backbone"
          else:
            target_selection = f"{target} and chain {target_chain} and resi {first_base_target_chain}-{last_base_target_chain} and backbone"
          #print("Target selection")
          #print(target_selection)
          target_residues = list(set([at.resi for at in cmd.get_model(target_selection).atom]))
          probe_residues = list(set([at.resi for at in cmd.get_model(moving_selection).atom]))

          #print("Check for missing atoms")
          #print(target_residues)
          #print(probe_residues)
          #Hanlde missing atoms

          #Check for missing atoms
          exclude_from_target, exclude_from_probe = [], []
          if len(target_residues) != len(probe_residues):
              print("Error: Different number of residues in target and probe selections")
          else:
              for o, target_res in enumerate(target_residues):
                  atom_names_target, atom_names_probe = [], []
                  for at in cmd.get_model(target_selection).atom:
                      if target_res == at.resi:
                          if not at.name in atom_names_target:
                            atom_names_target.append(at.name)
                  print("Atom names target")
                  print(atom_names_target)
                  for at in cmd.get_model(moving_selection).atom:
                      if probe_residues[o] == at.resi:
                          if not at.name in atom_names_probe:
                            atom_names_probe.append(at.name)
                  print("Atom names probe")
                  print(atom_names_probe)
                  for at_name_target in atom_names_target:
                      if not at_name_target in atom_names_probe:
                          print(f"Warning: {at_name_target} of {probe_residues[o]} not found in probe selection. Omitting from target selection")
                          exclude_from_target.append((target_res, at_name_target))
                  for at_name_probe in atom_names_probe:
                      if not at_name_probe in atom_names_target:
                          print(f"Warning: {at_name_probe} of {target_res} not found in target selection. Omitting from probe selection")
                          exclude_from_probe.append((probe_residues[o], at_name_probe))

              for u, item in enumerate(exclude_from_target):
                  if u == 0:
                      target_selection = f"{target_selection} and not ((i. {item[0]} and name {item[1]})"
                  elif u == len(exclude_from_target) - 1:
                      target_selection = f"{target_selection} or (i. {item[0]} and name {item[1]}))"
                  else:
                      target_selection = f"{target_selection} or (i. {item[0]} and name {item[1]})"

              for u, item in enumerate(exclude_from_probe):
                  if u == 0:
                      moving_selection = f"{moving_selection} and not ((i. {item[0]} and name {item[1]})"
                  elif u == len(exclude_from_probe) - 1:
                      moving_selection = f"{moving_selection} or (i. {item[0]} and name {item[1]}))"
                  else:
                      moving_selection = f"{moving_selection} or (i. {item[0]} and name {item[1]})"

          print(f"Aligning {moving_selection} on {target_selection}")
          cmd.select("moving", moving_selection)
          cmd.select("target", target_selection)
          cmd.color("red", moving_selection)
          cmd.color("blue", target_selection)

          try:
              align_result = cmd.pair_fit("moving", "target")
          except:
              align_result = None
          
          if align_result:
            clashes_protein = "clashes_protein"
            clashes_target_DNA = "clashes_DNA"
            clashes_protein_byres = "clashes_protein"
            clashes_target_DNA_byres = "clashes_DNA"
            # Calculate clashes between probe protein and target protein
            cmd.select(clashes_protein, f"({probe} and polymer.protein) within {cutoff} of ({target} and polymer.protein)")
            # Calculate clashes between probe protein and target DNA
            cmd.select(clashes_target_DNA, f"({probe} and polymer.protein) within {cutoff} of ({target} and polymer.nucleic)")
            clashes_protein_count = cmd.count_atoms(clashes_protein)
            print(f"Clashes with target protein (atom count): {clashes_protein_count}")
            clashes_DNA_count = cmd.count_atoms(clashes_target_DNA)
            print(f"Clashes with target DNA (atom count): {clashes_DNA_count}")
            # Count residue clashes
            cmd.select(clashes_protein_byres, f"{clashes_protein} byres {clashes_protein}")
            cmd.select(clashes_target_DNA_byres, f"{clashes_target_DNA} byres {clashes_target_DNA}")
            clashes_protein_res_count = cmd.count_atoms(f"{clashes_protein_byres} and name CA")
            clashes_DNA_res_count = cmd.count_atoms(f"{clashes_target_DNA_byres} and name CA")
            print(f"Clashes with target protein (residue count): {clashes_protein_res_count}")
            print(f"Clashes with target DNA (residue count): {clashes_DNA_res_count}")
            # Delete objects to avoid memory overflow
            cmd.delete(clashes_protein)
            cmd.delete(clashes_target_DNA)
            cmd.delete(clashes_protein_byres)
            cmd.delete(clashes_target_DNA_byres)

            if not align_both_strands:
              first_base_target_chain = str(first_base_target_chain).replace('\\', '')
              last_base_target_chain = str(last_base_target_chain).replace('\\', '')
            else:
              first_base_target_chain_1 = str(first_base_target_chain_1).replace('\\', '')
              last_base_target_chain_1 = str(last_base_target_chain_1).replace('\\', '')
              first_base_target_chain_2 = str(first_base_target_chain_2).replace('\\', '')
              last_base_target_chain_2 = str(last_base_target_chain_2).replace('\\', '')
            if align_both_strands:
              pose_name = f"{target_chain_1}_{first_base_target_chain_1}-{last_base_target_chain_1}_{target_chain_2}_{first_base_target_chain_2}-{last_base_target_chain_2}"
            else:
               pose_name = f"{target_chain}_{first_base_target_chain}-{last_base_target_chain}"

            # Write out models (nuclesome + current superimposed probe) for Rosetta scoring
            if write_models:
              model_name = f"{target}_{probe}_{pose_name}.pdb"
              editing.copy_to(model_name, f"({target} or {probe})")
              cmd.save(os.path.join(dir_name, model_name), model_name)
              cmd.delete(model_name)

            # Retain models with no. residue clashes less than or equal to 'clash_keep'
            if clashes_protein_res_count < int(clash_keep) + 1:
                editing.copy_to(f"probe_{pose_name}", probe)
                to_remove = f'probe_{pose_name} and not (({moving_selection.replace(f"{probe} and ", "")}) or polymer.protein)'
                print(to_remove)
                cmd.remove(to_remove)
                cmd.group(f"clash_keep_{clash_keep}_probe_length_{probe_length})", f"nonclashing + probe_{pose_name}")

            if align_both_strands:
              target_chain = f"{target_chain_1}_{target_chain_2}"
              #For logging only the forward strand is given
              first_base_target_chain = first_base_target_chain_1
              last_base_target_chain = last_base_target_chain_1

            print(f"DNA chain {target_chain} bases {first_base_target_chain}-{last_base_target_chain}: {probe} clashes with {clashes_protein_count} non-hydrogen atoms in {clashes_protein_res_count} residues")
            row = {'chain': target_chain,
                      'start': str(first_base_target_chain).replace('\\', ''),
                      'end': str(last_base_target_chain).replace('\\', ''),
                      'Atom clashes Protein/Protein': clashes_protein_count,
                      'Residue clashes Protein/Protein': clashes_protein_res_count,
                      'Atom clashes Protein/DNA': clashes_DNA_count,
                      'Residue clashes Protein/DNA': clashes_DNA_res_count,
                      'Residue clashes overall': clashes_protein_res_count + clashes_DNA_res_count,
                      'Atom clashes overall': clashes_protein_count + clashes_DNA_count,
                      'RMSD': align_result}
            results.append(row)

          # Catch superpostion failures
          else:
            print(f"DNA chain {target_chain_1} bases {first_base_target_chain} - {last_base_target_chain}: {probe} superpositon failed and clashes were not calculated"
                    f" – check for missing atoms, alternate conformations, or unsual bases at these positions")
            row = {'chain': target_chain,
                    'start': first_base_target_chain,
                      'end': last_base_target_chain,
                        'Atom clashes Protein/Protein': "-",
                    'Residue clashes Protein/Protein': "-",
                      'Atom clashes Protein/DNA': "-",
                        'Residue clashes Protein/DNA': "-",
                                  'RMSD': "-"}
            results.append(row)

      df = pd.DataFrame(results)
      print(df)
      csv_out = f'{target}_{probe}_probe_length_{probe_length}_clashes.csv'
      df.to_csv(csv_out)

    if only_plot:
        if os.path.exists(csv_out):
            df = pd.read_csv(csv_out)
        else:
           raise OSError("CSV not found!")
    df = df.dropna()
    df['Residue clashes Protein/Protein'] = pd.to_numeric(df['Residue clashes Protein/Protein'], errors='coerce')
    df['Residue clashes Protein/DNA'] = pd.to_numeric(df['Residue clashes Protein/DNA'], errors='coerce')
    df['Overall residue clashes'] = df['Residue clashes Protein/Protein'] + df['Residue clashes Protein/DNA']
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['RMSD'] = pd.to_numeric(df['RMSD'], errors='coerce')
    grouped_df = df.groupby('chain')
    for group_name, group_data in grouped_df:
      g, ax = plt.subplots()

      # Plot Residue Clashes data on the left y-axis
      group_data.plot(x='start', y='Residue clashes Protein/Protein', ax=ax, label='Protein/Protein')
      group_data.plot(x='start', y='Residue clashes Protein/DNA', ax=ax, label='Protein/DNA')
      group_data.plot(x='start', y='Overall residue clashes', ax=ax, label='Overall residue clashes')
      ax.set_xlabel('Start')
      ax.set_ylabel('Residue Clashes')
      ax.set_title(f'Residue Clashes for {group_name}')
      ax.legend(loc='upper left')

      ax2 = ax.twinx()

      ax2.plot(group_data['start'], group_data['RMSD'], color='red', label='RMSD', linestyle='--')
      ax2.set_ylabel('RMSD')
      ax2.legend(loc='upper right')

      plt.savefig(f'{target}_{probe}_residue_clashes_combined_{group_name}_probe_length_{probe_length}.png', format='png')
  except KeyboardInterrupt:
     sys.exit()
  except Exception as e:
     print(e)
     print(traceback.print_exc())
