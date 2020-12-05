#!/usr/bin/env python3

"""

Extract bonds from .mol file for building an atoms_list file for Blacomcalc

Author: (c) Matteo Paolieri, University of Cologne, Dec 2020
Version: 1.1
License: MIT

Docs: https://github.com/mtplr/blacomcalc

"""

import argparse
import numbers as num
import os


def change_btype(bt):
    # converts bond order to letter
    # 1 = s, 2 = d, 3 = t

    try:
        bt = int(bt)
    except Exception as e:
        print(f'Error with bond type conversion occurred: {e}')

    if bt == 1:
        bt = 's'
    elif bt == 2:
        bt = 'd'
    elif bt == 3:
        bt = 't'
    elif bt == 4:
        print('Found a bond order 4 ("aromaticity"). Not supported. Please check your bond types.')
    else:
        print('Problem with the bond type conversion.')

    return bt


def extract(molfile, n):

    with open(molfile, 'r') as f:
        data = f.readlines()  # read .mol file just to see first where is the bonds part

    i = 0
    start_line = 0

    for line in data:

        line = line.split()

        if len(line) == 7:  # the turning point of the mol file is when it goes to 7 entries = bonds stuff
            start_line = i  # at what line start to read the file
            break
        else:
            i += 1
            continue  # if the line doesn't have 7 entries, skip it

    with open(molfile, 'r') as f:  # this time at the right point
        data = f.readlines()[start_line:]  # read all the data from the mol file

    del data[-1]  # remove the last item containing "END", in .mol file

    with open(f'bonds-{molfile[:-4]}.txt', 'w') as f:  # remove ".mol" ext

        # need to read the line every 3 chars because the limit is 999 atoms for .mol files

        print(f'#M', file=f)

        i = 0

        for line in data:

            i += 1  # keep trace of line number

            chars = []  # for every line reset the chars

            for ch in line:
                if ch == ' ':
                    ch = ''
                    chars.append(ch)
                elif ch.isdigit() is True:  # write spaces and integers only
                    chars.append(ch)
                else:
                    continue

            a1 = chars[0] + chars[1] + chars[2]
            a2 = chars[3] + chars[4] + chars[5]
            bt = chars[8]

            new_bt = change_btype(bt)

            print(f'{a1} {a2} {new_bt}', file=f)

            if (i % n) == 0:  # check if the number of molecules is a multiple of the line number; add #M delimiter
                if i == len(data):  # check if we are at the last '#M'
                    print(f'#M', file=f)
                else:
                    print(f'#M\n#M', file=f)
            else:
                continue


def valid_file(param):  # check if the file is .mol
    base, ext = os.path.splitext(param)
    if ext.lower() not in '.mol':
        print(f'Error: file must have a .mol extension. Used {param} instead.')
    return param


def main (mol, n):
    try:
        extract(mol, n)
    except Exception as e:
        print(f'Error: {e}')


if __name__ == "__main__":
    # parser for shell
    parser = argparse.ArgumentParser(description='Extract bond information from .mol file and make it '
                                                 'readable for Blacomcalc atoms_list file.')
    parser.add_argument('mol_file', type=valid_file, help="Input .mol file")
    parser.add_argument('number_atoms', type=int, help="Input: number of atoms of each molecule, e.g. "
                                                       "for 5 molecules of H2O: n = 3. The file will be splitted in "
                                                       "5 parts delimited by #M.")
    args = parser.parse_args()
    main(args.mol_file, args.number_atoms)
