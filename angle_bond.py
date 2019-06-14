#!/usr/bin/python
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from molmod.constants import boltzmann
from molmod.io.xyz import XYZFile
from molmod.ic import bend_angle
from molmod.ic import bond_length
from scipy.optimize import curve_fit
import json


def main(file_name, start_step, end_step, temp):
    """
    Loads geometries, generates file with time evolution
    of given angles and bonds. Prints force parameters for a harmonic oscillator
    that reproduces the behavior of the each bond/angle
    """

    # Store the trajectory as a numpy array
    xyz_file = XYZFile(file_name)
    geometries = xyz_file.geometries[start_step:end_step]

    # Read atom list from input file (input is 1-based indexing)  
    with open('atoms.txt', 'r') as f:
        atoms_input = json.load(f)
    atoms = [list(np.array(a) - 1) for a in atoms_input]

    # Calculate bonds and angles
    time = np.arange(geometries.shape[0]) + 1
    
    bonds = [get_bonds(geometries, i) for i in atoms if len(i) == 2]
    angles = [get_angles(geometries, i) for i in atoms if len(i) == 3]

    # Calculate histograms for bonds and angles and save them
    if bonds:
        for i, colvar in enumerate(bonds):
                all_distr, coefficients = generate_histogram(colvar, temp)
                name = convertLabel(atoms_input[i - 1])
                np.savetxt(name + ".hist", all_distr)
                np.savetxt(name + ".dat", colvar)
                np.savetxt(name + ".coeff", coefficients)
        #np.savetxt("all_bonds.dat", np.vstack((time, np.stack(bonds))).transpose())
        labels = [convertLabel(i) for i in atoms_input if len(i) == 2]
        all_bonds = pd.DataFrame(data=np.stack(bonds).transpose(), index=time, columns=labels)
        all_bonds.to_csv("all_bonds.dat", sep='\t')

    if angles:
        for i, colvar in enumerate(angles):
                all_distr, coefficients = generate_histogram(colvar, temp)
                name = convertLabel(atoms_input[i - 1])
                np.savetxt(name + ".hist", all_distr)
                np.savetxt(name + ".dat", colvar)
                np.savetxt(name + ".coeff", coefficients)
        #np.savetxt("all_angles.dat", np.vstack((time, np.stack(angles))).transpose())
        labels = [convertLabel(i) for i in atoms_input if len(i) == 3]
        all_bonds = pd.DataFrame(data=np.stack(angles).transpose(), index=time, columns=labels)
        all_bonds.to_csv("all_bonds.dat", sep='\t')


def generate_histogram(colvar, temp):
    # Define the histogram and shift the bins to have data on the centres
    hist, bin_edges = np.histogram(colvar, bins=50, density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    histogram = np.stack((bin_centres, hist))
    gaussian, probdistr, coefficients = fit_gaussian(histogram, temp)
    all_distr = np.stack((bin_centres, hist, gaussian, probdistr), axis=1)
    return all_distr, coefficients


def convertLabel(colvar):
    '''
    generates a label with the list of atoms  [70, 170] --> Bond_70_170
    '''
    if len(colvar) == 3:
        label = "Angle_"
    else:
        label = "Bond_"
    lab = '_'.join(str(x) for x in colvar)
    return label + lab


def get_bonds(geometries, atoms):
    '''
    This functions takes an array with the geometries, and a list with a pair of atoms 
    and returns the bonds evolution during the simulation
    '''
    number_of_steps = geometries.shape[0]
    number_of_atoms = geometries.shape[1]
    distance = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        distance[frame] = bond_length(geometries[frame, atoms])[0]
    return distance


def get_angles(geometries, atoms):
    '''
    This functions takes an array with the geometries, and a list with a triplet of atoms 
    and returns the angle evolution during the simulation

    '''
    number_of_steps = geometries.shape[0]
    number_of_atoms = geometries.shape[1]
    angle = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        angle[frame] = bend_angle(geometries[frame, atoms])[0]
    return angle


def fit_gaussian(data, temp, p0=[2, 0.1]):
    '''
    Takes histogram and bins of same shape and returns a fitted gaussian distribution
    and the force constants of the corresponding harmonic oscillator 
    '''
    bins, hist = data[0], data[1]
    coeff, var_matrix = curve_fit(gaussian_distribution, bins, hist, p0=p0)
    
    # Fit a gaussian distribution to the bond/angle distribution. 
    # We fit this distribution and not directly the oscillator distribution because sometimes this fit fails.
    gauss_fit = gaussian_distribution(bins, *coeff)

    # Obtain force constant for the oscillator that generates this distribution
    #append it to the coefficient list. 
    # sigma = sqrt(k/2pi kb T)
    kb = boltzmann  # kb in hartree
    k = (temp * kb) / (coeff[1]**2)
    all_coefficients = np.append(coeff, k)
    
    # Check if k replicates the real distribution
    coeff_distr = k, coeff[0], temp
    distr_fit = oscillator_distribution(bins, *coeff_distr)
    return gauss_fit, distr_fit, all_coefficients


def gaussian_distribution(x, *p):
    '''
    Define a gaussian function to fit data. A = 1 / np.sqrt(2 * np.pi * sigma**2)
    '''
    mu, sigma = p
    return (1 / (np.sqrt(2 * np.pi) * sigma)) * \
        np.exp(-(x - mu)**2 / (2. * sigma**2))


def oscillator_distribution(x, *p):
    '''
    Define the probability distribution function at given k and T for testing purposes.
    This distribution should be the same as the gaussian one.
    '''
    k, x0, temp = p
    kb = boltzmann
    return(1 / (np.sqrt(2 * np.pi * kb * temp / k)) * \
        np.exp(-(k / (2 * kb * temp)) * (x - x0)**2))

def plot_everything(hist, colvar, coefficients):





if __name__ == "__main__":
    msg = "angle_bond -p <path/to/trajectory> -st <start frame> -et <end frame> -t <temperature>"
    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-p', required=True, help='path to the xyz trajectory')
    parser.add_argument(
        '-st',
        required=False,
        default=0,
        type=int,
        help='starting time of the simulation (default=0)')
    parser.add_argument('-et', required=False, default=-
                        1, type=int, help='ending time of the simulation')
    parser.add_argument('-t', required=False, default=298, help='temperature')
    args = parser.parse_args()
    main(args.p, args.st, args.et, args.t)
