#!/usr/bin/python
import numpy as np
import molmod as md
import argparse
from functools import partial
from molmod.io.xyz import XYZFile
from molmod.ic import bend_angle
from molmod.ic import bond_length
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy.stats import norm

def main(file_name)#, start_step, end_step):
    xyz_file = XYZFile(file_name)
    outputName = "bond_angles.dat"
    outFilenameGraph = "bond_anglesGraph"
    frames = xyz_file.geometries
    atomsi = [[70,170],[165,68],[170,70,9],[170,70,167],[165,68,5],[165,68,267]]  #add this in the command line
    atoms = [ list(np.array(a)-1) for a in atomsi ]

    bonds  = map(partial(get_bonds, frames), filter(lambda at: len(at) == 2, atoms))
    angles = map(partial(get_angles, frames), filter(lambda at: len(at) == 3, atoms))
    numbers = [np.arange(frames.shape[0])+1]

    allThing = np.concatenate((numbers, bonds, angles))
    np.savetxt(outputName,np.transpose(allThing))
    plotBonds_Angles(outFilenameGraph,allThing,atomsi)

def plotBonds_Angles(outFilename,allThing,atomsi):
    ''' 
    Save files with data as a function of time and histogram
    '''
    (num,steps) = allThing.shape
    for i in range(num)[1:]:
        # Define the histogram and shift the bins to have data on the centres
        hist, bin_edges = np.histogram(allThing[i], bins=50,density=True)
        bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
        histogram = np.stack((bin_centres, hist))
        name = convertLabel(atomsi[i-1])
        gaussian, coefficients = fit_gaussian(histogram)
#        mu, std = norm.fit(allThing[i])
#        gaussian2 = norm.pdf(bin_centres, mu, std)
        printall = np.stack((bin_centres, hist, gaussian), axis =1 )
        np.savetxt(name+".hist", printall)
        np.savetxt(name+".dat", allThing[i])
        np.savetxt(name+".coeff", coefficients)

def convertLabel(xs):
    ''' 
    generates a label with the list of atoms  [70, 170] --> Bond_70_170
    '''
    if len(xs) == 3:
       label = "Angle_"
    else:
       label = "Bond_"
    lab = '_'.join(str(x) for x in xs)
    return label + lab

def get_bonds(frames,atoms): 
    '''
    This functions wants an array with the frames, and a list with a pair of atoms
    '''
    number_of_steps = frames.shape[0]
    number_of_atoms = frames.shape[1]
    distance = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        distance[frame] = bond_length(frames[frame,atoms])[0]*0.529177
    return distance

def get_angles(frames,atoms): 
    '''
    This functions wants an array with the frames, and a list with a triplet of atoms
    '''
    number_of_steps = frames.shape[0]
    number_of_atoms = frames.shape[1]
    angle = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        angle[frame] = bend_angle(frames[frame,atoms])[0]#*(180./np.pi)
    return angle

def fit_gaussian(data, p0=[2,0.1]):
    '''
    takes histogram and bins of same shape
    '''
    bins, hist = data[0],data[1] 
    coeff, var_matrix = curve_fit(gauss, bins, hist, p0=p0)
    # Get the fitted curve
    hist_fit = gauss(bins, *coeff)
    # get k of the oscillator that generates this distribution and append it to the coefficient list
    # sqrt(k/2pi kb T)
    temp = 351
    kb = 1.38064852E-23
    JtoHartree = 2.293710449e17
    kb_Hartree = 3.1668114E-6
    k = (temp*kb_Hartree)/coeff[1]
#    k = (temp*kb*JtoHartree)/coeff[1]
    allcoeff = np.append(coeff,k)
#    print('k = {0} Hartree/Bohr**2,  r0 = {1} Bohr'.format(k, mu))
    return hist_fit, allcoeff

def gauss(x, *p):
    '''
    Define a gaussian function to fit data. A = 1/np.sqrt(2*np.pi*sigma**2)
    '''
    mu, sigma = p
    return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-(x-mu)**2/(2.*sigma**2))
    #return (1/np.sqrt(2*np.pi*sigma**2))*numpy.exp(-(x-mu)**2/(2.*sigma**2))
    
if __name__ == "__main__":
    msg = " angle_bond -p <path/to/trajectory>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-p', required=True, help='path to the xyz trajectory')
    parser.add_argument('-st', required=False, help='starting time of the simulation')
    parser.add_argument('-et', required=False, help='ending time of the simulation')
#    parser.add_argument('-a', required=False, help='give three values for an angle')
#    parser.add_argument('-c', required=False, help='give two values for a bond length')
    args = parser.parse_args()
    main(args.p)
#    main(args.p, args.st, args.et)
