# Blacomcalc

A simple Computational Chemistry Python script to calculate bond lengths,
BLA value (Bond Length Alternation) and center of mass (CoM) of the given 
`atoms_list` molecules and bonds, starting from an `.xyz` standard file.

**Features:** 

* Calculate every bond length distance between the inserted atoms, for the given molecule(n)

* Calculate the BLA according to the definition: the used formula 
is the average of single bond distances minus the average of double bonds

* Calculate the center of mass of the desired inserted molecules

* Calculate the distance between the desired center of masses (**WIP**)

# Usage

## To use it

Clone this repo:

`git clone https://github.com/mtplr/blacomcalc` 

Launch it:

```bash
calc_bla.py molecule.xyz atoms_list > output-blacomcalc
```

## input file

The script reads a standard `.xyz` file containing the molecule geometry and an input file `atoms_list` containing:

* The first row: assign one number to the total number of molecules

* `#M` delimiter, then two atoms for each row, the ones you want the program to calculate the distance between, 
and a third column that specifies the bond type, everything separated by a space. This can be repeated with 
multiple molecules separated by `#M`. For the definition of the single and double bonds the user will decide (before) 
which is single and which is double based on the chemical structure given by the **neutral form**. 

* At the end of molecule blocks `#COM` must be inserted, followed by the couple of molecules you want to 
calculate the distance between their center of masses.

```
2
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
```

The example above means: "calculate the distance between atoms 14 and 13, which corresponds to 
a double bond" and so on, for two molecules. Finally, calculate the center of mass between the molecule 1 and 2.

`s` indicates a single bond and `d` a double.

## .xyz file

It can be generated with Avogadro, Molden, Gaussian...

Example of file generated with Avogadro:

```
25

C         -2.03043       -7.02776        1.91131
C         -3.02001       -6.11648        1.55490
C         -2.55205       -4.85898        1.17223
C         -0.74691       -6.48160        1.80939
S         -0.96570       -4.92699        1.27622
C         -3.36434       -3.78577        0.76351
...       ...             ...            ...
```

where the number `25` followed by an empty line indicates the number of atoms in the molecule, followed by an empty line.

## Example

If you launch the example files here enclosed (`pyrl.xyz` and `atoms_list`) it is possible to compare 
the calculated bond lengths obtained here (`outputbla`) against e.g. the Avogadro output 
[(picture here)](BL_avogadro.png).

# License

Author: (c) Matteo Paolieri, University of Cologne

Version: 30.11.2020

License: MIT (see LICENSE)