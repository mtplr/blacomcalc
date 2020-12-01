#!/usr/bin/env python3

"""

Blacomcalc

A simple Computational Chemistry Python script to calculate bond lengths,
BLA value (Bond Length Alternation) and center of mass (CoM) of the given
atoms_list molecules and bonds, starting from an .xyz standard file.

Author: (c) Matteo Paolieri, University of Cologne
Version: 30.11.2020
License: MIT

"""


import argparse
import numpy as np
import math
import re as regex
import time


def calc_distance(x1, x2, y1, y2, z1, z2):
    # calculate the distance of two 3D points (x,y,z)

    dist: float = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return dist


def define_mass(m):
    # assign a mass value to the specified atom
    # (source IUPAC https://www.ciaaw.org)

    masses_dic = {
        "H": 1.0078250322,
        "C": 12.0,
        "N": 14.003074004,
        "O": 15.994914619,
        "S": 31.972071174,
        "P": 30.973761998,
        "F": 18.998403163,
        "Cl": 34.9688527
    }

    return masses_dic[m]


def calc_com(xyz_file, atom_min, atom_max):

    # calculate the Center of Mass (com)
    # vR = (1/M)sum(mi*vr_i) where v = vector

    # the algorithm search max and min value in the atoms_list file,
    # and retrieve all atoms in the xyz
    # this way it is possible to know what is what
    # and the range for each molecule

    atom_min: int = atom_min-1  # remove 1 because of the 0 index in lists
    atom_max: int = atom_max

    with open(xyz_file) as xyz:

        # read the xyz file and skip the first two lines
        xyz = xyz.readlines()[2:]
        coord = xyz[atom_min:atom_max]

    # generate a vector with vectors
    # containing all atoms, e.g.
    # [['C', '-2.03043', '-7.02776', '1.91131'], [], ...]

    list_of_atoms = []

    tot_atoms = range((len(coord)))

    for line in tot_atoms:
        atom_i = coord[line].split()
        list_of_atoms.append(atom_i)

    sum_mi_ri = np.array([0.0, 0.0, 0.0])
    sum_mi = 0

    for atom in list_of_atoms:

        # assign a mass to the atom type
        # then calculate the sum of masses
        label = atom[0]
        mi = float(define_mass(label))
        sum_mi = sum_mi + mi

        # remove label, so I have (x,y,z) vector
        atom.pop(0)

        # change every value to float and put it in a numpy vector
        vr_i = np.array(atom, dtype=np.float32)

        # dot product mass*r_i
        mi_ri = np.dot(mi, vr_i)

        # sum all over m_i*r_i
        sum_mi_ri = sum_mi_ri + mi_ri

    # get the coordinates of CoM
    # (and covert np.array to list)
    com = (np.dot(1/sum_mi, sum_mi_ri)).tolist()

    x = com[0]
    y = com[1]
    z = com[2]

    com_coord = [x,y,z]

    return com_coord


def calc_bla(xyz_file, bla_data):
    # open .xyz file and put all columns in a vector
    # then calculate ALL desired bond lengths and print them
    # then calculate the BLA value

    # initialize atom vectors

    label_list = []
    x_list = []
    y_list = []
    z_list = []

    # generate vectors with labels and coord

    with open(xyz_file) as coord:

        # skip the first 0,1 lines of xyz file (start from 2)
        coord = coord.readlines()[2:]

        for line in range(len(coord)):

            # put every line of xyz file in a separated vector
            column = coord[line].split()

            # take 'line' of column 'n'
            label = column[0]
            x = column[1]
            y = column[2]
            z = column[3]

            # append 'line' of column 'n' into a vector
            label_list.append(label)
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)

    # split the string into lines (list)
    # and read bla data from the second line
    atoms_to_search = bla_data.split('\n')

    # init index vectors
    atom1_list = []
    atom2_list = []
    bond_type_list = []

    for atoms in atoms_to_search:

        # obtain data in atoms_list
        # to each index it corresponds a bond type
        # i.e. s (single), d (double), t (triple)

        column = atoms.split()
        atom1 = column[0]
        atom2 = column[1]
        bond_type = column[2]

        atom1_list.append(atom1)
        atom2_list.append(atom2)
        bond_type_list.append(bond_type)

    # number of atoms in atoms_list
    # for each couple there's an index number linked with the vector
    # this is used to retrieve the bond length info later
    tot_atoms = range(len(atom1_list))

    # init vector with all bond lengths
    bls_list = []

    for index in tot_atoms:
        # generate atom index from bla file

        # retrieve atom number
        # remove 1 because it starts from 0
        atom_num1 = int(atom1_list[index])-1
        atom_num2 = int(atom2_list[index])-1

        # retrieve all info from vectors and calc distance
        # index is the same for each couple of atoms

        a1 = label_list[atom_num1]
        a2 = label_list[atom_num2]
        x1 = float(x_list[atom_num1])
        x2 = float(x_list[atom_num2])
        y1 = float(y_list[atom_num1])
        y2 = float(y_list[atom_num2])
        z1 = float(z_list[atom_num1])
        z2 = float(z_list[atom_num2])
        bond_type = bond_type_list[index]

        bond_length = calc_distance(x1, x2, y1, y2, z1, z2)
        bls_list.append(bond_length)

        # note: add again +1 because the 0 is not counted here
        print(f'Calculated bond length between atoms '
              f'{str(a1)} ({str(atom_num1+1)}) and {str(a2)} ({str(atom_num2+1)}), {bond_type}: \n '
              f'{str(bond_length)}')

    # finally, calculate the BLA value:
    # dBLA = sum(d_single)/N – sum(d_double)/M
    # using bond_type_list vector, I retrieve all bond lengths in bls_list
    # to calculate with the above formula

    s = 0          # sum of all single bond lengths
    s_count = 0    # number of single bonds
    d = 0          # sum of all double bond lengths
    d_count = 0    # number of double bonds

    for i in range(len(bond_type_list)):

        bond_type = bond_type_list[i]   # get the bond type (s,d)
        sum_bonds = float(bls_list[i])  # get the bond length

        if bond_type == 's':
            s = s + sum_bonds
            s_count += 1
        elif bond_type == 'd':
            d = d + sum_bonds
            d_count += 1

    bla_value = s / s_count - d / d_count

    print(f'\n\n#####################\n'
          f'Total single bonds: {s_count}\n'
          f'Total double bonds: {d_count}\n'
          f'Final BLA value is (Å): {bla_value}\n'
          f'#####################\n\n')

    # write a BLA.dat file to plot it with gnuplot
    with open("BLA.dat", "w") as f:  # create a .dat file (! it overwrites !)
        print(f'# bond\tbond length (Å)', file=f)
        n = 1
        for line in bls_list:
            print(f'{n}\t{line}', file=f)  # print all lengths + bond index
            n += 1


def parse(input_file):
    # split and parse the atoms_list file

    with open(input_file) as input_file:

        text = input_file.read()

        try:
            n_of_molecules = int(text[0])
        except:
            print("Error in the atoms_list file (number of molecules)")

        # find all molecule blocks (max split = number of molecules)
        molecule_blocks = text.split('\n#M\n')
        # find the block related to CoM calculation
        com_blocks = text.split('#COM\n')

        # cleaning
        molecule_blocks.pop(0)  # remove first blank space
        molecule_blocks.pop(len(molecule_blocks)-1) # remove last entry with COM
        com_blocks.pop(0)  # remove first blank space
        com_blocks = com_blocks[0].split(' ')   # just take the first entry and split it, otherwise is e.g. ['1 2']

    # create a tuple for the output
    parsed_text = (n_of_molecules, molecule_blocks, com_blocks)

    return parsed_text


def find_min_max_in_str(string):
    # find max and min given a string

    numbers = regex.findall(r'\d+', string)  # \d+ = find all numbers

    # convert all str in int in numbers
    for i in range(len(numbers)):
        numbers[i] = int(numbers[i])

    max_n = max(numbers)
    min_n = min(numbers)

    results = (min_n, max_n)

    return results


def main(xyz, input_file):

    start_time = time.time()

    # first, parse the text

    parsed_text = parse(input_file)

    n_of_molecules = int(parsed_text[0])
    molecule_blocks = parsed_text[1]
    com_block = parsed_text[2]  # what molecule(s) for CoM calc

    # find for what atoms and molecules calculate the CoM
    # for every molecule in com_block I search the corresponding
    # min and max atom index in the corresponding atoms_list block
    # so that I can pass it to calc_com function and select
    # the right atoms in the xyz file

    com_atoms_min = []
    com_atoms_max = []

    for molecule in com_block:

        molecule_string = (molecule_blocks[int(molecule) - 1])

        min_atom = find_min_max_in_str(molecule_string)[0]
        max_atom = find_min_max_in_str(molecule_string)[1]

        com_atoms_min.append(min_atom)
        com_atoms_max.append(max_atom)

    # for each molecule...
    for molecule_bla in range(n_of_molecules):

        print(f'---------------------\n\n'
              f'BLA for MOLECULE: {molecule_bla+1}  \n\n'
              f'---------------------\n\n')

        # calculate bond lengths and BLA value
        calc_bla(xyz, molecule_blocks[molecule_bla])

    coms_coord_list = []  # initialize

    for molecule_com in range(len(com_block)):

        # calculate Center of Mass of the single molecules,
        # specifying the range of every molecule
        com_coord = calc_com(xyz, int(com_atoms_min[molecule_com]), int(com_atoms_max[molecule_com]))

        x = com_coord[0]
        y = com_coord[1]
        z = com_coord[2]

        coms_coord_list.append(com_coord)  # list with all the CoM's coordinates

        print(f'-------------------------------------------------\n'
              f'The center of mass of the molecule no. {int(molecule_com)+1} is at:\n'
              f'x = {x}\n'
              f'y = {y}\n'
              f'z = {z}\n'
              f'-------------------------------------------------\n\n')

    # calculate the distance between the selected CoM's

    # TODO calculate the distance between multiple molecules !!

    # print time
    print(f'######## Completed in: {round(((time.time() - start_time)*1000),2)} ms #########')


if __name__ == "__main__":
    # parser for shell
    parser = argparse.ArgumentParser(description='A small Python script to calculate BLA value and bond distances.',
                                     epilog='Usage: blacomcalc.py xyz_file bla_file. Output for data is BLA.dat')
    parser.add_argument('xyz_file', type=str, help="Input .xyz file: first 2 rows must be skipped!")
    parser.add_argument('input_file', type=str, help="Input .bla file containing the desired distances between "
                                                   "atoms. It suffices to write two atoms and the type of bond "
                                                   "(s = single, d = double) in columns, "
                                                   "all separated by a space, e.g.: \"atom1 atom2 d\". "
                                                   "First line must begin with a comment "
                                                   "(although # is not mandatory).")
    args = parser.parse_args()
    main(args.xyz_file, args.input_file)