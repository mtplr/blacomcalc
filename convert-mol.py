#!/usr/bin/env python3

"""

Extract bonds from .mol file for building an atoms_list file for Blacomcalc

Author: (c) Matteo Paolieri, University of Cologne, Dec 2020
Version: 1.0
License: MIT

"""

import argparse


def change_btype(bt):
    # converts bond order to letter
    # 1 = s, 2 = d, 3 = t

    try:
        bt = int(bt)
    except:
        print("Error with bond type conversion occurred.")

    if bt == 1:
        bt = 's'
    elif bt == 2:
        bt = 'd'
    elif bt == 3:
        bt = 't'
    elif bt == 4:
        print('Found a bond order 4. Not supported.')
    else:
        print('Problem with the bond type conversion.')

    return bt


def extract(molfile):

    with open(molfile, 'r') as f:
        data = f.readlines()

    with open('bonds-from-mol.txt', 'w') as f:

        for line in data:

            column = line.split()

            if len(column) == 7:  # the turning point of the mol file is when it goes to 7 entries = bonds stuff
                a1 = column[0]
                a2 = column[1]
                bt = column[2]

                new_bt = change_btype(bt)

                print(f'{a1} {a2} {new_bt}', file=f)

            else:
                continue  # if the line doesn't have 7 entries, skip it


def main (mol):
    extract(mol)


if __name__ == "__main__":
    # parser for shell
    parser = argparse.ArgumentParser(description='Extract bond information from .mol file and make it '
                                                 'readable for Blacomcalc atoms_list file.')
    parser.add_argument('mol_file', type=str, help="Input .mol file")
    args = parser.parse_args()
    main(args.mol_file)
