#!/usr/bin/python
import numpy as np
from molmod.io.xyz import XYZFile
import argparse
import os


def main(file1_name, file2_name, start_step, end_step, trim, unitcell_name):
    """
    Takes a xyz trajectory and can combine it with the ptcv trajectory, trim it
    and add unit cell parameters in the title lines so that PLUMED can process it
    if the unit cell is variable

    """
    if os.path.isfile(file1_name):
        xyz_file1 = XYZFile(file1_name)
        xyz_file3 = xyz_file1
        geo1 = xyz_file1.geometries[start_step:end_step]

        if file2_name is not None:
            if os.path.isfile(file2_name):
                xyz_file2 = XYZFile(file2_name)
                geo2 = xyz_file2.geometries[start_step:end_step]

                xyz_file2.symbols = ('He',)
                min_size = min(geo1.shape[0], geo2.shape[0])

                xyz_file3.geometries = np.concatenate(
                    (geo1[:min_size:trim], geo2[:min_size:trim]), axis=1)
                xyz_file3.symbols = np.concatenate(
                    (xyz_file1.symbols, xyz_file2.symbols), axis=0)
                xyz_file3.write_to_file('coord-ptcv.xyz')

                if unitcell_name is not None:
                    unitcellmerge(
                        xyz_file3,
                        start_step,
                        end_step,
                        trim,
                        unitcell_name)
            else:
                print('%s does not exist' % (file2_name))
        else:
            xyz_file3.geometries = geo1[::trim]
            if unitcell_name is not None:
                xyz_file3 = unitcellmerge(
                    xyz_file3, start_step, end_step, trim, unitcell_name)
            xyz_file3.write_to_file('coord-processed.xyz')
    else:
        print('%s does not exist' % (file1_name))


def unitcellmerge(xyz_file, start_step, end_step, trim, unitcell_name):
    """
    Takes a md-CELLFILE.cell in cp2k format and a xyz trajectory object from molmod
    returns a trajectory with cell parameters in the comment line for each step
    """
    with open(unitcell_name, 'r') as f:
        output = f.readlines()[start_step:end_step]
    cell = [x.split()[2:11] for x in output if not x.startswith('#')]
    cell2 = [' '.join(x) for x in cell[::trim]]
    xyz_file.titles = cell2
    return xyz_file


if __name__ == "__main__":
    msg = "merge-traj.py -i <path/to/trajectory> -tr <trim (print every n steps)> -st <start frame> -et <end frame> -uc <path/to/cp3k-unit-cell>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument(
        '-i',
        required=False,
        default='mof-pos-1.xyz',
        help='path to the xyz trajectory of file 1')
    parser.add_argument('-p', required=False, default=None,
                        help='path to the xyz trajectory of file 2')
    parser.add_argument(
        '-st',
        required=False,
        default=0,
        type=int,
        help='starting time of the simulation (default=0)')
    parser.add_argument('-et', required=False, default=-
                        1, type=int, help='ending time of the simulation')
    parser.add_argument(
        '-tr',
        required=False,
        default=1,
        type=int,
        help='trim every X steps')
    parser.add_argument(
        '-uc',
        required=False,
        default=None,
        help='unit cell file')
    args = parser.parse_args()
    main(args.i, args.p, args.st, args.et, args.tr, args.uc)
