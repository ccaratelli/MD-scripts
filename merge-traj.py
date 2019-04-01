#!/usr/bin/python
import numpy as np
from molmod.io.xyz import XYZFile
import argparse
import os

def main(file1_name, file2_name, start_step, end_step, trim):
    if os.path.isfile(file1_name):
        xyz_file1 = XYZFile(file1_name)
        xyz_file3 = XYZFile(file1_name)
        geo1= xyz_file1.geometries[start_step:end_step]

        if file2_name != None:
            if os.path.isfile(file2_name):
                xyz_file2 = XYZFile(file2_name)
                geo2 = xyz_file2.geometries[start_step:end_step]
        
                xyz_file2.symbols = ('He',)
                min_size = min(geo1.shape[0], geo2.shape[0])
        
                xyz_file3.geometries = np.concatenate((geo1[:min_size:trim], geo2[:min_size:trim]), axis=1)
                xyz_file3.symbols =  np.concatenate((xyz_file1.symbols, xyz_file2.symbols), axis=0)
                xyz_file3.write_to_file('coord-ptcv.xyz')
            else:
                print('%s does not exist' %(file2_name))
        else:
            xyz_file3.geometries = geo1[::trim]
            xyz_file3.write_to_file('coord-reduced.xyz')
    else:
        print('%s does not exist' %(file1_name))


if __name__ == "__main__":
    msg = " angle_bond -i <path/to/trajectory> -tr trim -st start_time -et end_time"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-i', required=False, default='mof-pos-1.xyz', help='path to the xyz trajectory of file 1')
    parser.add_argument('-p', required=False, default=None, help='path to the xyz trajectory of file 2')
    parser.add_argument('-st', required=False, default=0, type=int, help='starting time of the simulation (default=0)')
    parser.add_argument('-et', required=False, default=-1, type=int, help='ending time of the simulation')
    parser.add_argument('-tr', required=False, default=1, type=int, help='trim every X steps')
    args = parser.parse_args()
    main(args.i, args.p, args.st, args.et, args.tr)
