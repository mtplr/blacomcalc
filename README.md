[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](http://makeapullrequest.com)
[![GitHub issues](https://img.shields.io/github/issues/mtplr/blacomcalc)](https://github.com/mtplr/blacomcalc/issues)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

# Blacomcalc

A simple script written in Python to calculate bond lengths, BLA value (Bond Length Alternation), center of masses (CoM), distances between center of masses and bond angles of the given molecules, starting from a standard `.xyz` molecular geometry file.

In addition, a `BLA.csv` file containing all bond lengths is created for each molecule. It can be plotted with e.g. Gnuplot, OriginLab or Matplotlib.

Please, read carefully how to properly create the **input file**, [here](#input-file).

## Credits

Author: © **Matteo Paolieri**, University of Cologne, 2020.

Special thanks for critics and suggestions to: Dr. Daniele Fazzi, Nora Gildemeister, Ludovico Pavesi.

License: **MIT** (see LICENSE).

## Contributing

Suggestions and contributions are always welcome! 😉 Please discuss larger changes (e.g. adding torsion angles) via issue before submitting a pull request.

## Table of Contents

- [Blacomcalc](#blacomcalc)
  - [Credits](#credits)
  - [Contributing](#contributing)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Usage](#usage)
  - [Input file](#input-file)
  - [Get data from .mol files](#get-data-from-mol-files)
  - [.xyz file](#xyz-file)
  - [Quickstart and examples](#quickstart-and-examples)
    - [Example of input files](#example-of-input-files)
    - [Quickstart](#quickstart)

## Features

* Calculation of every bond length distance between the inserted atoms, for the given molecules

* Calculation of the BLA according to the definition: the used formula, as reported below, is the average of single bond distances minus the average of double bonds. It is also possible to average specific bonds defined by the user.

* Creation of a `BLA.csv` file for every molecule. This can be plotted with e.g. Gnuplot, OriginLab, Matplotlib

* Calculation of the center of mass (CoM) of the desired inserted molecules

* Calculation of the distance between the desired center of masses

* Calculation of the bond angle for the three specified atoms

* Extraction of what bonds are formed between what atoms, and the bond type (single, double, triple), from a `.mol` file using the converter script `convert-mol.py`. The output is ready for Blacomcalc `input_file` (see [later](#get-data-from-mol-files-wip)).
  
## Usage

Clone this repo:

```bash
git clone https://github.com/mtplr/blacomcalc
```

Launch the script:

```bash
blacomcalc.py molecule.xyz input_file.txt > output-blacomcalc.txt
```

Help command:

```bash
$  blacomcalc.py -h

usage: blacomcalc.py [-h] xyz_file input_file

Blacomcalc (v. 2.0.10). A simple script to calculate bond lengths, BLA values (Bond
Length Alternation), bond angles, center of masses (CoM) and distance between the CoMs
of the given molecules, starting from a standard .xyz molecular geometry file.

positional arguments:
  xyz_file    Geometry input file (.xyz). First 2 rows are skipped.
  input_file  Input file (plain text). See documentation.

optional arguments:
  -h, --help  show this help message and exit

(c) Matteo Paolieri 2020, License MIT. Docs: https://github.com/mtplr/blacomcalc
```

## Input file

The script reads a standard `.xyz` file (the extension is **mandatory**!) containing the molecule geometry and the indications on what to looking for are within an `input_file`. The extension of such a file doesn't matter, it is plain text, so no extension or e.g. `.txt`, `.inp` are equally fine.

Data must be reported according to the
following specifications:

* The first row: assign one number to the total number of molecules

* `#M` tag, then two atoms for each row, the ones you want the program to calculate the distance between, and a third column that specifies the bond type, everything separated by a space. This can be repeated with multiple molecules separated by `#M`. For the definition of the single and double bonds the user will decide (before) which is single and which is double based on the chemical structure given by the **neutral form**.

  * To quickly find the atoms involved in bonds, the [software Avogadro](https://avogadro.cc) can be used. It suffices to go on _Selection Tool (F11) > Selection Mode > Molecule_ and then click on settings on _Label_ in _Display Types_, and finally click to _Display only selective primitives_, as [shown here](img/avogadro_selection.png).

	* A new tool `convert-mol.py` to simplify this task has also been added 
	 (see [this section here](#get-data-from-mol-files-wip))

	* Using the above converter will give you only a list of bonds in the Blacomcalc format, but of course you can modify the input file by hand as you wish (bond types included)

* Every block must be specified within the appropriate tags: `#M...#M`, `#COM...#COM` or `#ANGLES...#ANGLES`

* Angles are calculated specifying between what atoms, e.g.: for `12 13 14` the angle centered on atom 13, which is
in the middle between 12 and 14 will be calculated
  
* At the end of the input file it is possible to leave a blank line or write `end` for completion.

* Bonds must be written in the input file as you expect the sorting in the final `BLA.csv` plot file, i.e.: if the single bond between atoms 1 and 2 `1 2 s` in the input is in position `14` of the molecule block, it will be the 14th bond in the final `BLA.csv` plot file.
  
* If needed, it is possible to average specific bonds. For example, due to a resonating double bond among near atoms. In that case, those will be counted as **one** single or double bond for BLA calculation. It suffices to write `avs` (that means *average bond calculation to sum up to single bonds*) followed by a number. They must be written **at the end of the molecule block** and **sequentially**, this way they will be replaced in the `BLA.csv` file with final averaged bonds. 

  For example, in context:

    ```text
    #M
    12 13 s
    11 12 avs1
    12 33 avs1 
    7 8 avs2
    45 55 avs2
    ```

  Here, bonds with the same `avs` number will be averaged (i.e. those defined by atoms 11, 12 and 12, 33 for `avs1` and those between 7, 8 and 45, 55 for `avs2`) and counted as single bonds for the final BLA calculation. Those bonds length will be printed in the output, but removed from the BLA.csv files. There, it will be reported only the averaged value for each `avs` number, **sequentially, at the end of file**.

  In that example above, bonds `11 12` and `12 33` are removed from `.csv` file and replaced with their average. The same for `avs2`.

  An identical function for this purpose but related to the summation of double bonds needs to be implemented (e.g. `avd`). Collaborations are welcome!

* `COM` and `ANGLES` blocks can be empty, just specify the tags between the word `null`. Be careful, because it is **case-sensitive**!

  For example, if you don't want to calculate the center of mass and/or the angles, just write the blocks as:

  ```text
  #COM
  null
  #COM
  ```

  ```text
  #ANGLES
  null
  #ANGLES
  ```

* **Caveat!** Center of mass (CoM) is calculated **only with respect to the provided bonds**! For instance, if you insert only the bonds related to the molecule's backbone, Blacomcalc will calculate only the center of mass related **only to those atoms!**. It can be a good approximation (for example to calculate distance between molecules), but if you want the "full" CoM, you have to insert all the molecules bonds in the `input_file`.

## Get data from .mol files

A quick tool to extract bonds information in the way Blacomcalc reads it, is implemented: just use **convert-mol.py**. This is just a raw script, so it must still be considered **WIP**!

This feature is currently work in progress, and it is reports data for Blacomcalc only for **oligomers** of the same molecule. However, it can be used for different molecules only if `#M` delimiters are put later by hand.

You can launch it as:

```bash
convert-mol.py file.mol atoms-number
```

Where `atoms-number` is the number of atoms of **one** molecule of the oligomer. Put it as e.g. 999 to not put delimiters (e.g. with three molecules in an oligomer: `convert-mol.py geometry.mol 3`).

It may be useful to first convert a `.xyz` file in a `.mol` file using eg. Avogadro or [OpenBabel](http://openbabel.org/wiki/Main_Page) first, and then run the script.

This way it is possible to extract useful information to build the `input_file` file, since it contains the couple of atoms for each bond and it specifies what kind of bond (bond order) it is involved (1 = single, 2 = double, 3 = triple, 4 = aromacity).

The script converts it in `1 = s, 2 = d, 3 = t`.

`.mol` files can present also a number 4 for bond type, i.e. _aromaticity_. This is not supported with `convert-mol.py` thus is needed to change them if necessary. A check for this in the script is present.

The `.mol` file possesses this turning point:

```text
...
   -4.9566   -2.6804    1.3230 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9603   -3.5558    2.9163 H   0  0  0  0  0  0  0  0  0  0  0  0
   -5.1669   -2.3024   -0.3827 S   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3744   -3.5300    3.2906 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 11  1  0  0  0  0
  2 10  1  0  0  0  0
  3  6  1  0  0  0  0
...
```

Here, e.g. `1  2  1  0  0  0  0` means: "single bond formed by atoms 1 and 2", and so on. The script parses this in order to get such lines in a form `1 2 s`.

The complete anatomy of a `.mol` file can be found [here, on Wikipedia](https://en.wikipedia.org/wiki/Chemical_table_file#Molfile).

## .xyz file

It can be generated with Avogadro, Molden, Gaussian etc.

Example of `.xyz` file (generated with Avogadro):

```text
25

C         -2.03043       -7.02776        1.91131
C         -3.02001       -6.11648        1.55490
C         -2.55205       -4.85898        1.17223
C         -0.74691       -6.48160        1.80939
S         -0.96570       -4.92699        1.27622
C         -3.36434       -3.78577        0.76351
...       ...             ...            ...
```

where the number `25` indicates the number of atoms in the molecule, followed by a line.

## Quickstart and examples

### Example of input files

In the [quickstart folder](quickstart) there are two examples of input files and their application.
Here there two "abstract" examples.

```text
4
#M
14 13 d
15 16 d
12 13 s
12 15 s
12 10 s
10 11 s
. . .
#M
24 25 d
35 36 s
52 53 s
. . .
#M
#COM
1 2
1 4
3 2
. . .
#COM
#ANGLES
12 13 14
199 1 53
. . .
#ANGLES
end
```

The example above means: *calculate the distance between atoms 14 and 13, which corresponds to a double bond* and so on, for two molecules. Finally, calculate the center of mass between the molecule 1 and 2 and the bond angle between the atoms 12 13 and 14.

`s` indicates a single bond and `d` a double.

This other one, instead, will calculate only the BLA value for one molecule, and one bond angle between three atoms.

```text
1
#M
14 13 s
15 16 d
12 13 s
12 15 s
12 10 d
10 11 s
14 199 avs1
188 5 avs1
. . .
#M
#COM
null
#COM
#ANGLES
12 13 14
38 1 187
#ANGLES
end
```

### Quickstart

```bash
blacomcalc.py molecule.xyz input_file.txt > blacomcalc-output.txt
```

If you launch the example files here enclosed in the _Quickstart_ folder, with `molecule.xyz` as geometry and `input_file.txt` as the input file, it is possible for example to compare the calculated bond lengths obtained here against e.g. the Avogadro molecule output [(picture here)](img/BL_avogadro.png).

[GNF2-xTB](https://github.com/grimme-lab/xtb) was used for molecular optimizations.

The overall printed outputs are reported in `output-blacomcalc.txt`.
