#!/usr/bin/env python3

"""

    +============================================================================+
    |                                                                            |
    |   (  _`\ (_ )                                               (_ )           |
    |   | (_) ) | |    _ _    ___    _     ___ ___     ___    _ _  | |    ___    |
    |   |  _ <' | |  /'_` ) /'___) /'_`\ /' _ ` _ `\ /'___) /'_` ) | |  /'___)   |
    |   | (_) ) | | ( (_| |( (___ ( (_) )| ( ) ( ) |( (___ ( (_| | | | ( (___    |
    |   (____/'(___)`\__,_)`\____)`\___/'(_) (_) (_)`\____)`\__,_)(___)`\____)   |
    |                                                                            |
    |   Matteo Paolieri, University of Cologne, 2020                             |
    +============================================================================+


A simple Computational Chemistry Python script to calculate bond lengths,
BLA value (Bond Length Alternation) center of mass (CoM), distance between
the center of masses, and bond angles of the given molecules and bonds,
starting from an .xyz standard file.

Author: (c) Matteo Paolieri, University of Cologne, Dec 2020
Version: 2.0.5
License: MIT

Docs: https://github.com/mtplr/blacomcalc

"""


import argparse
import os
import collections, functools, operator
import numpy as np
import math
import re as regex
import time


def calc_distance(x1, x2, y1, y2, z1, z2):
    # calculate the distance of two 3D points (x,y,z)

    dist: float = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return dist


def find_min_max_in_str(string):
    # find max and min given a string

    numbers = regex.findall(r'\d+', string)  # \d+ = find one or more numbers

    # convert all str in int
    for i in range(len(numbers)):
        numbers[i] = int(numbers[i])

    max_n = max(numbers)
    min_n = min(numbers)

    results = (min_n, max_n)

    return results


def define_mass(m):
    # assign a mass value to the specified atom
    # (source IUPAC https://www.ciaaw.org)

    masses = {
        "H": 1.0078250322,
        "C": 12.0,
        "N": 14.003074004,
        "O": 15.994914619,
        "S": 31.972071174,
        "P": 30.973761998,
        "F": 18.998403163,
        "Cl": 34.9688527
    }

    return masses[m]


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

    com_coord = [x, y, z]

    return com_coord


def calc_bla(xyz_file, bla_data, molecule_number):
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

        column = atoms.split(" ")
        atom1 = column[0]
        atom2 = column[1]
        this_bond_type = column[2]

        atom1_list.append(atom1)
        atom2_list.append(atom2)
        bond_type_list.append(this_bond_type)

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
        this_bond_type = bond_type_list[index]

        bond_length = calc_distance(x1, x2, y1, y2, z1, z2)
        bls_list.append(bond_length)

        # note: add again +1 because the 0 is not counted here
        print(f'Calculated bond length between atoms '
              f'{str(a1)} ({str(atom_num1+1)}) and {str(a2)} ({str(atom_num2+1)}), {this_bond_type}: \n '
              f'{str(bond_length)}')

    # finally, calculate the BLA value:
    # dBLA = sum(d_single)/N – sum(d_double)/M
    # using bond_type_list vector, I retrieve all bond lengths in bls_list
    # to calculate with the above formula

    s = 0          # sum of all single bond lengths
    s_count = 0    # number of single bonds
    d = 0          # sum of all double bond lengths
    d_count = 0    # number of double bonds

    avs_count = 0  # 'avs' = average of bonds to be added to single bonds count for BLA calculation
    avs_list_key = []
    avs_list_bond_length = []

    for i in range(len(bond_type_list)):  # It includes the average between selected bonds

        this_bond_type = bond_type_list[i]   # get the bond type (s,d, t, avs_n...)
        this_bond_length = float(bls_list[i])  # get the bond length

        if this_bond_type == 's':
            s = s + this_bond_length
            s_count += 1
        elif this_bond_type == 'd':
            d = d + this_bond_length
            d_count += 1
        elif this_bond_type.find('avs') >= 0:  # Thanks to Ludovico Pavesi for the useful comments

            # TODO Do this also for avd, i.e. for averaging bonds that need to be summed within double bonds

            # split this_bond_type with respect to the relative number
            # and make an average between the selected atoms
            # and consider them as single bonds

            avs_list_key.append(this_bond_type[3:])  # remove 'avs' and put as a key to average the same ones
            avs_list_bond_length.append(this_bond_length)

    # first average all the bonds you need to
    # and add it to the single bonds sum

    if len(avs_list_key) > 0:  # check if there are any bond to be average before calculating the BLA

        avs_dic_list = []  # init new LIST of dictionaries

        for i in range(0, len(avs_list_key)):

            dic = {avs_list_key[i]:avs_list_bond_length[i]}  # generate a key-value dictionary
            avs_dic_list.append(dic)  # append to the dictionary: now we have a set of key/values and need to avg them

        avs_key_list_no_repetitions = list(dict.fromkeys(avs_list_key))  # remove duplicates in avs_list_key

        # calculate average for each key (i.e. every avs1 is averaged
        # then every avs2, etc... All of these are counted in single bonds)

        for key in list(dict.fromkeys(avs_list_key)):
            avs_key_average = sum(item.get(str(key), 0) for item in avs_dic_list) / len(avs_key_list_no_repetitions)
            s_count += 1
            s = s + avs_key_average
            avs_count += 1

    # calculate the final BLA value:
    bla_value = s/s_count - d/d_count

    print(f'\n=========================================================\n'
          f'Molecule no.: {molecule_number}\n'
          f'Total single bonds: {s_count}\n'
          f'Total double bonds: {d_count}\n'
          f'Total averaged bonds, counted as single: {avs_count}\n'
          f'Final BLA value is (Å): {bla_value}\n'
          f'=========================================================\n')

    # write a BLA.dat file to plot it with gnuplot
    with open(f'BLA-{molecule_number}.dat', "w") as f:
        print(f'# bond no. bond length (Å)', file=f)
        n = 1
        for line in bls_list:
            print(f'{n} {line}', file=f)  # print all lengths + bond index
            n += 1


def parse(input_file):
    # split and parse the atoms_list file

    with open(input_file) as input_file:
        text = input_file.read()

    # find the number (int) of molecules from the first line
    n_of_molecules = int(text[0])

    # split the blocks
    molecule_blocks = text.split('\n#M\n')
    com_blocks = text.split('\n#COM\n')
    angle_blocks = text.split('\n#ANGLES\n')

    # generate list of molecule blocks
    molecule_blocks = molecule_blocks[1:]  # remove no. of molecules
    molecule_blocks.pop(n_of_molecules)  # remove last block (stuff after the last molecule)

    # generate a list of string of CoM block
    # e.g. [['1', '2'], ['3', '4']]
    com_blocks = com_blocks[1]
    com_blocks = com_blocks.split('\n')

    parsed_com_blocks = []

    for molecule in com_blocks:
        j = molecule.split(' ')
        parsed_com_blocks.append(j)

    com_blocks = parsed_com_blocks

    # generate a list of strings of atoms, for angles block, to calculate the angle between
    # e.g. [['14', '13', '15'], ['16', '17', '18']]

    angle_blocks = angle_blocks[1]
    angle_blocks = angle_blocks.split('\n')

    parsed_angle_blocks = []

    for molecule in angle_blocks:
        j = molecule.split(' ')
        parsed_angle_blocks.append(j)

    angle_blocks = parsed_angle_blocks

    # create a tuple for the output
    # (type: int, lst, str, lst)
    parsed_text = (n_of_molecules, molecule_blocks, com_blocks, angle_blocks)

    return parsed_text


def calc_distance_coms(com_block, coms_coord):

    # calculate the distance between the selected CoM's SEQUENTIALLY

    # com_block = vector with the number of molecules you want to
    # calculate the CoM between, eg ['1', '2', '1', '2']

    # coms_coord = vector with the coordinates of each CoM
    # e.g. [[x1,y1,z1], [x2,y2,z2], ...]

    # first convert all str to int
    for i in range(0, len(com_block)):
        com_block[i] = int(com_block[i])

    try:
        n_of_molecules = int(len(com_block))
        if n_of_molecules == 1:
            print('Distance between center of masses was not calculated. There\'s only one molecule.\n')

        elif (n_of_molecules % 2) == 0 or 1:  # if the number is even or there's only one couple...

            # ...calc distance sequentially e.g. "10 12 13 14", first 10-12 then 13-14
            # for loop, iterating by 2
            for i in range(0, n_of_molecules, 2):
                mol1 = i
                mol2 = i+1

                coord1 = coms_coord[mol1]
                coord2 = coms_coord[mol2]

                x1 = coord1[0]
                y1 = coord1[1]
                z1 = coord1[2]

                x2 = coord2[0]
                y2 = coord2[1]
                z2 = coord2[2]

                dist = calc_distance(x1, x2, y1, y2, z1, z2)

                print(f'The distance between the center of masses of molecules {int(com_block[mol1])} and'
                      f' {int(com_block[mol2])} is:\n{dist} Å\n\n')

    except Exception as e:
        print(f'ERROR in calculating distances between molecules. Please, control #COM input.\n{e}')


def calc_bond_angles(xyz_file, angle_block):

    # first convert all str to int, so we have e.g. [1, 2, 3, 4, 5, 6]
    # -1 must be done because I need the index starting from 0,
    # in order to read the xyz file and retrieve all the coordinates
    # the goal is to calculate SEQUENTIALLY the angle between 1, 2, 3, then 4, 5, 6...

    for i in range(0, len(angle_block)):
        angle_block[i] = int(angle_block[i])-1

    # find the coordinates

    with open(xyz_file) as xyz:

        # read the xyz file and skip the first two lines
        xyz = xyz.readlines()[2:]

    # generate a vector with vectors
    # containing all atoms, e.g.
    # [['C', '-2.03043', '-7.02776', '1.91131'], [], ...]

    list_of_atoms = []

    for line in angle_block:
        atom_i = xyz[line].split()
        list_of_atoms.append(atom_i)

    # build a coordinate vector for each atom

    coord_atoms = []
    labels = []

    for atom in list_of_atoms:
        label = atom[0]  # TODO: add labels in printing
        x = atom[1]
        y = atom[2]
        z = atom[3]

        coord = (x, y, z)

        labels.append(label)
        coord_atoms.append(coord)

    # given the coordinates, calculate the bond angle
    # in pseudo-code: angle = arccos[dot_product(v1, v2)/((norm(v1)*norm(v2))]

    j = 0  # init to count the iterations

    for i in range(0, int(len(angle_block)), 3):

        j += 1  # to count the iterations

        # find indices
        atom1 = i
        atom2 = i+1
        atom3 = i+2

        # find labels
        label1 = labels[atom1]
        label2 = labels[atom2]
        label3 = labels[atom3]

        # find vectors v1(a,b) v2(b,c)
        a = np.array(coord_atoms[atom1], dtype=float)
        b = np.array(coord_atoms[atom2], dtype=float)
        c = np.array(coord_atoms[atom3], dtype=float)

        v1 = np.array([a[0] - b[0], a[1] - b[1], a[2] - b[2]], dtype=float)
        v2 = np.array([c[0] - b[0], c[1] - b[1], c[2] - b[2]], dtype=float)

        # calculate the angle in degrees

        v1v2 = np.vdot(v1, v2)
        nv1 = np.linalg.norm(v1)
        nv2 = np.linalg.norm(v2)
        angle = round(math.degrees(math.acos(v1v2/(nv1*nv2))), 2)

        print(f'\n'
              f'The angle no. {j} between {label1} {label2} {label3} is: {angle}°'
              f'\n')


def main(xyz, input_file):

    # TODO: clean a little bit the spaghetti-code here and put everything in the right function

    print('''\n
 +============================+
 |                            |
 |   Blacomcalc OUTPUT file   |
 |   v. 2.0.5                 |
 |                            |
 +============================+
    \n\n''')

    try:
        start_time = time.time()

        # first, parse the text

        parsed_text = parse(input_file)

        n_of_molecules = int(parsed_text[0])  # number of molecules
        molecule_blocks = parsed_text[1]  # what molecules for BLA

        # ============= BLA ================================================

        # for each molecule...
        for molecule_bla in range(n_of_molecules):

            molecule_number_bla = int(molecule_bla) + 1

            print(f'\n\n---------------------\n\n'
                  f'BLA for MOLECULE: {molecule_number_bla}  \n\n'
                  f'---------------------\n\n')

            # calculate bond lengths and BLA value
            calc_bla(xyz, molecule_blocks[molecule_bla], molecule_number_bla)

        # ============= COM ================================================

        com_block = parsed_text[2]  # what molecules for CoM (lst of couples)

        if 'null' in com_block[0]:

            print("No COM to calculate.")

        else:

            # find for what atoms and molecules calculate the CoM
            # for every molecule in com_block I search the corresponding
            # min and max atom index in the corresponding atoms_list block
            # so that I can pass it to calc_com function and select
            # the right atoms in the xyz file

            # search each molecule in the COM block and find max and min value
            # to get the atoms range in the xyz file

            com_atoms_min = []
            com_atoms_max = []

            # first, join the lst of lst (trick by Alex Martelli, 2009)
            molecules_for_com = [item for sublist in com_block for item in sublist]

            for single_molecule_com in molecules_for_com:
                # I need to pass a list of molecules (com_block) as [1,2,3...n]

                # then calculate min and max for each molecule

                molecule_string = (molecule_blocks[int(single_molecule_com) - 1])

                min_atom = find_min_max_in_str(molecule_string)[0]
                max_atom = find_min_max_in_str(molecule_string)[1]

                com_atoms_min.append(min_atom)
                com_atoms_max.append(max_atom)

            # now search the selected molecule ranges and get the needed atoms!

            coms_coord_list = []  # initialize

            for molecule_com in range(len(molecules_for_com)):

                # calculate the Center of Mass of the single molecules,
                # specifying the range of every molecule
                com_coord = calc_com(xyz, int(com_atoms_min[molecule_com]), int(com_atoms_max[molecule_com]))

                x = com_coord[0]
                y = com_coord[1]
                z = com_coord[2]

                coms_coord_list.append(com_coord)  # list with all the CoM's coordinates

                molecule_number = int(molecules_for_com[molecule_com])

                print(f'\n'
                      f'The center of mass of the molecule no. {molecule_number} is at:\n'
                      f'x = {x}\n'
                      f'y = {y}\n'
                      f'z = {z}\n'
                      f'\n')

            # TODO: print the CoM's coordinate of the selected molecules without repetitions

            # TODO: add torsional angle calculation

            # calculate the distance between the selected CoM's
            # using all the molecules sequentially
            calc_distance_coms(molecules_for_com, coms_coord_list)

        # ============= ANGLES =============================================

        angle_block = parsed_text[3]  # what molecules for angles (lst of str, triples)

        if 'null' in angle_block[0]:

            print("No bond angle to calculate.")

        else:

            # first join the lst of lst
            atoms_for_angles = [item for sublist in angle_block for item in sublist]

            # calculate bond angle
            if atoms_for_angles[0] == '0':  # dirty fix to check empty vector. TODO: add a clean fix
                print(f'No angles to calculate. {atoms_for_angles[0]}')
            else:
                print(f'\n====================\n'
                      f'ANGLES'
                      f'\n====================\n')
                calc_bond_angles(xyz, atoms_for_angles)

        # ============= TIMING ===============================================

        # print time
        print(f'\n######## Completed in: {round(((time.time() - start_time)*1000), 2)} ms #########')

    except IndexError as e:

        print(f'ERROR: {e}.\nThere might be a problem with the number of atoms in the input file.'
              f' Maybe with line 1 or CoM\'s in atoms_list')

    except Exception as e:

        print(f'ERROR: {e}.')  # TODO enhance error handling...


def valid_file(param):  # check if the file is .xyz
    base, ext = os.path.splitext(param)
    if ext.lower() not in '.xyz':
        print(f'ERROR: file must have a .xyz extension. Used {param} instead.')
    return param


if __name__ == "__main__":
    # parser for shell
    parser = argparse.ArgumentParser(description='A simple Computational Chemistry Python script to '
                                                 'calculate bond lengths, BLA value (Bond Length Alternation), '
                                                 'bond angles and center of mass (CoM) '
                                                 'of the given atoms_list molecules and bonds, '
                                                 'starting from a standard .xyz file.',
                                     epilog='Usage: blacomcalc.py xyz_file atoms_list. Output for plotting '
                                            'data is BLA-n.dat, with n = number of molecule in atoms_list.')
    parser.add_argument('xyz_file', type=valid_file, help="Input .xyz file: first 2 rows must be skipped!")
    parser.add_argument('input_file', type=str, help="Input atoms_list file containing the desired distances between "
                                                     "atoms. See README.")
    args = parser.parse_args()
    main(args.xyz_file, args.input_file)
