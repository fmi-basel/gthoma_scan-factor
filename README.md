# Scan Factor

## Description

A PyMOL plugin to iteratively align the DNA of a probe object (e.g. TF bound to DNA motif) onto the DNA of a target object (e.g. nucleosome) and calculate the clashes between protein atoms of probe and DNA and/or protein atoms of target.

## Installation

Requires PyMOL with python > 3.x

In the PyMOL console install dependencies
`pip install pandas`
`pip install matplotlib`

Download the scan_factor.py file and install it as a plugin as described under (https://pymolwiki.org/index.php/Plugins).

## Usage example

```
fetch 6t90, async=0
create NCP, 6t90 and c. A+B+C+D+E+F+G+H+I+J
create OCTSOX, 6t90 and (c. K+L or (c. I and i. 3-24) or (c. J and i. 124-145))
scan_factor NCP, OCTSOX, probe_center=True
```

## Contributors

Richard D. Bunker and Georg Kempf