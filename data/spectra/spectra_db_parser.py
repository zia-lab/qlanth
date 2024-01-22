#!/usr/bin/env python3

import os
import json
import numpy as np

elst_fname = 'elst.dbf'
spin_fname = 'spin.dbf'
crystf_fname = 'recr.dbf'
export = True

def parse_spectra_dbf():
    counter = 0
    matrix_elements = []
    op_index_to_name = {1: 'F2',
    2: 'F4',
    3: 'F6',
    4: 'α',
    5: 'β',
    6: 'γ',
    7: 't2',
    8: 't3',
    9: 't4',
    10: 't6',
    11: 't7',
    12: 't8'
    }
    with open(elst_fname,'rb') as file:
        while True:
            bytes_to_read = [4,6,2,4][counter%4]
            data = file.read(bytes_to_read)
            if data == b'':
                break
            if counter % 4 == 0:
                header = data
                n = header[0]
                op_index = op_index_to_name[header[2]]
            elif counter % 4 == 1:
                # try:
                #     deg = int(matrix_element[2])
                #     print(deg)
                # except:
                #     extra = file.read(1)
                #     print(extra)
                matrix_element = data.decode('utf-8')
                braterm, keterm = matrix_element[:3], matrix_element[3:]
                braterm =  braterm.strip()
                keterm = keterm.strip()
            elif counter % 4 == 2:
                data = ''
            elif counter % 4 == 3:
                matrix_value = float(np.frombuffer(data, dtype=np.float32)[0])
            # print(data)
            if counter % 4 == 3:
                matrix_elements.append(((braterm, keterm), (n, op_index), matrix_value))
            counter += 1 

    electro_matrix_elements = matrix_elements

    counter = 0
    matrix_elements = []
    elstops = {1: 'ζ',
    2: 'M0',
    3: 'M2',
    4: 'M4',
    5: 'P2',
    6: 'P4',
    7: 'P6'}
    with open(spin_fname,'rb') as file:
        while True:
            bytes_to_read = [5,7,4][counter%3]
            data = file.read(bytes_to_read)
            if data == b'':
                break
            if counter % 3 == 0:
                header = data
                n = header[0]
                twoJay = header[2]
                op_index = header[-1]
                op = elstops[op_index]
            elif counter % 3 == 1:
                matrix_element = data[:-1].decode('utf-8')
                braterm, keterm = matrix_element[:3], matrix_element[3:]
                braterm =  braterm.strip()
                keterm = keterm.strip()
            elif counter % 3 == 2:
                matrix_value = float(np.frombuffer(data, dtype=np.float32)[0])
            if counter % 3 == 2:
                matrix_elements.append(((braterm, keterm), (n, twoJay, op), matrix_value))
            counter += 1 

    magnetic_matrix_elements = matrix_elements

    counter = 0
    matrix_elements = []
    with open(crystf_fname,'rb') as file:
        while True:
            bytes_to_read = [8,4,4,4][counter%4]
            data = file.read(bytes_to_read)
            if data == b'':
                break
            if counter % 4 == 0:
                matrix_element = data[1:-1].decode('utf-8')
                braterm, keterm = matrix_element[:3], matrix_element[3:]
                braterm =  braterm.strip()
                keterm = keterm.strip()
            elif counter % 4 == 1:
                U2val = float(np.frombuffer(data, dtype=np.float32)[0])
            elif counter % 4 == 2:
                U4val = float(np.frombuffer(data, dtype=np.float32)[0])
            elif counter % 4 == 3:
                U6val = float(np.frombuffer(data, dtype=np.float32)[0])
            # print(data)
            if counter % 4 == 3:
                matrix_elements.append(((braterm, keterm), U2val, U4val, U6val))
            counter += 1 

    crystal_field_uk_matrix_elements = matrix_elements

    return magnetic_matrix_elements, crystal_field_uk_matrix_elements, electro_matrix_elements

if __name__ == '__main__':
    magnetic_matrix_elements, crystal_field_uk_matrix_elements, electro_matrix_elements = parse_spectra_dbf()
    if export:
        print("exporting")
        with open('spectra_matrix_elements.py', 'w') as file:
            file.write('''
# This file contains the matrix elements used in Spectra, an old code from Argonne National Lab.

# The matrix elements are stored in the following variables:
# magnetic_matrix_elements in the form of ((braTerm, ketTerm), (num_electrons, twoJay, operator_label), matrixElement)
# crystal_field_uk_matrix_elements in the form of ((braTerm, ketTerm), U2val, U4val, U6val)
# electro_matrix_elements in the form of ((braTerm, ketTerm), (num_electrons, operator_label), matrixElement)

# In all cases a missing element is implied to be zero.\n\n''')
            file.write('magnetic_matrix_elements = [')
            for line in magnetic_matrix_elements:
                file.write(str(line) + ',\n')
            file.write(']\n\n')
            file.write('crystal_field_uk_matrix_elements = [')
            for line in crystal_field_uk_matrix_elements:
                file.write(str(line) + ',\n')
            file.write(']\n\n')
            file.write('electro_matrix_elements = [')
            for line in electro_matrix_elements:
                file.write(str(line) + ',\n')
            file.write(']\n\n')
