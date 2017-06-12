import numpy as np
import molmod as md
import argparse
#from molmod.io.xyz import XYZReader
from functools import partial
from molmod.io.xyz import XYZFile
from molmod.ic import bend_angle
from molmod.ic import bond_length
import matplotlib as mpl
import matplotlib.pyplot as plt
def main(file_name):
#    xyz = XYZReader(file_name)
#    mols = [[m for m in f ] for f in xyz ]
#    frames = [j for i,j in mols]
#    tensor = np.asarray(frames)
#    number_of_steps = tensor.shape[0]
#    number_of_atoms = tensor.shape[1]
# # initialization of array
#    bonds = np.empty(number_of_steps,number_of_atoms,number_of_atoms)
#    angles = np.empty(number_of_steps)    
#    for i in range(number_of_steps): 
#        for j in range(number_of_atoms):
#        bond[i]  = bond_length(tensor[i,[atom1,atom2])[0]
#        angle[i] = bend_angle(tensor[i,[atom1,atom2,atom3])[0]

    xyz_file = XYZFile(file_name)
    outputName = "bond_angles.dat"
    outFilenameGraph = "bond_anglesGraph"
    frames = xyz_file.geometries
#    atomsi = [[170,70,9],[170,70,167],[70,170]]  #add this in the command line
    atomsi = [[165,68,5],[165,68,267],[165,68]]  #add this in the command line
    atoms = [ list(np.array(a)-1) for a in atomsi ]

    bonds  = map(partial(get_bonds, frames), filter(lambda at: len(at) == 2, atoms))
    angles = map(partial(get_angles, frames), filter(lambda at: len(at) == 3, atoms))
    numbers = [np.arange(frames.shape[0])+1]
    allThing = np.concatenate((numbers, bonds, angles))
    np.savetxt(outputName,np.transpose(allThing))
    plotBonds_Angles(outFilenameGraph,allThing,atomsi)
#    list_of_quantitities = filter(lambda at: len(at) == 2 and len(at) == 3, atoms)
    
#    for i in atoms:
#        if len(atoms[i]) == 2:
#           bond = get_bonds(frames,atoms[i])
#            list_of_quantities.append(bond)
#        if len(atoms[i]) == 3:
#            angle = get_angles(frames,atoms[i])
#            list_of_quantities.append(angle)
def plotBonds_Angles(outFilename,allThing,atomsi):
    (num,steps) = allThing.shape
    for i in range(num)[1:]:
        tupleZ = zip(*np.histogram(allThing[i], bins=50))
        name = convertLabel(atomsi[i-1])
        np.savetxt(name, tupleZ)

def convertLabel(xs):
    if len(xs) == 3:
       label = "Angle_"
    else:
       label = "Bond_"
    lab = '_'.join(str(x) for x in xs)
    return label+lab+".hist"

def get_bonds(frames,atoms): 
    """
    This functions wants an array with the frames, and a list with a pair of atoms
    """
    number_of_steps = frames.shape[0]
    number_of_atoms = frames.shape[1]
    distance = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        distance[frame] = bond_length(frames[frame,atoms])[0]
    return distance

def get_angles(frames,atoms): 
    """
    This functions wants an array with the frames, and a list with a triple of atoms
    """
    number_of_steps = frames.shape[0]
    number_of_atoms = frames.shape[1]
    angle = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        angle[frame] = bend_angle(frames[frame,atoms])[0]
    return angle

    
if __name__ == "__main__":
    msg = " angle_bond -p <path/to/trajectory>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-p', required=True, help='path to the xyz trajectory')
#    parser.add_argument('-a', required=False, help='give three values for an angle')
#    parser.add_argument('-c', required=False, help='give two values for a bond length')
    args = parser.parse_args()
    main(args.p)
