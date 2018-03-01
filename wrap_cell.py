#!/usr/bin/python
from molmod.io.xyz import XYZFile
from molmod.units import angstrom
from molmod.unit_cells import UnitCell
import argparse
import numpy as np

def main(file_name):
    xyz_file = XYZFile(file_name)
    frames = xyz_file.geometries
#    titles = xyz_file.titles

    # Unit cell, decide how to make this general 
    matrix=np.array([\
         [1.4731497044857509E+01, 3.2189795740722255E-02, 4.5577626559295564E-02],\
         [4.2775481701113616E-02, 2.1087874593411915E+01, -2.8531114198383896E-02],\
         [6.4054385616337750E-02, 1.3315840416191497E-02, 1.4683043045316882E+01]])
    matrix *= angstrom
    cell = UnitCell(matrix) 
    frac = UnitCell.to_fractional(cell, frames)

    xmin = np.min(frac[:,:,0])
    ymin = np.min(frac[:,:,1])
    zmin = np.min(frac[:,:,2])
    frac[:,:,0] -= -0.5#xmin
    frac[:,:,1] -= -0.5#ymin
    frac[:,:,2] -= -0.5#zmin
    decimals = np.modf(frac)[0] 
#    decimals[:,:,0] += xmin
#    decimals[:,:,1] += ymin
#    decimals[:,:,2] += zmin
    frac_wrapped = np.where(decimals < 0, 1 + decimals, decimals) 
#    frac_wrapped[:,:,0] += xmin
#    frac_wrapped[:,:,1] += ymin
#    frac_wrapped[:,:,2] += zmin
    cart_wrapped = UnitCell.to_cartesian(cell, frac_wrapped)

    xyz_file.geometries = cart_wrapped
    xyz_file.write_to_file(file_name.rsplit(".",1)[0]+"_wrapped.xyz")

if __name__ == "__main__":
    msg = " wrap_cell -p <path/to/trajectory>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-p', required=True, help='path to the xyz trajectory')
    args = parser.parse_args()
    main(args.p)
