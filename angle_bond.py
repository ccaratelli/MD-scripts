#!/usr/bin/python
import numpy as np
import argparse
import matplotlib.pyplot as plt
from molmod.constants import boltzmann
from molmod.units import angstrom
from molmod.io.xyz import XYZFile
from molmod.ic import bend_angle
from molmod.ic import bond_length
from scipy.optimize import curve_fit
import json


def main(file_name, start_step, end_step):
    """
    Loads geometries, generates file with time evolution
    of given angles and bonds.
    """
    # Store the trajectory as a numpy array
    xyz_file = XYZFile(file_name)
    geometries = xyz_file.geometries[start_step:end_step]

    # Read atom list from input file
    with open('atoms.txt', 'r') as f:
        atoms_input = json.load(f)
    atoms = [list(np.array(a) - 1) for a in atoms_input]

    # Calculate bonds and angles
    time = np.arange(geometries.shape[0]) + 1
    bonds = [get_bonds(geometries, i) for i in atoms if len(i) == 2]
    angles = [get_angles(geometries, i) for i in atoms if len(i) == 3]

    # Store all bonds and angles in an array
    all_variables = np.vstack((time, np.stack(bonds), np.stack(angles)))
    np.savetxt("bond_angles.dat", np.transpose(all_variables))

    n_colvars = all_variables.shape[0]
    for i in range(1, n_colvars):
        all_distr, coefficients = generate_histogram(all_variables[i])
        name = convertLabel(atoms_input[i - 1])
        np.savetxt(name + ".hist", all_distr)
        np.savetxt(name + ".dat", all_variables[i] / angstrom)
        np.savetxt(name + ".coeff", coefficients)


def generate_histogram(colvar):
    # Define the histogram and shift the bins to have data on the centres
    hist, bin_edges = np.histogram(colvar, bins=50, density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    histogram = np.stack((bin_centres, hist))
    gaussian, probdistr, coefficients = fit_gaussian(histogram)
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
    This functions wants an array with the geometries, and a list with a pair of atoms
    '''
    number_of_steps = geometries.shape[0]
    number_of_atoms = geometries.shape[1]
    distance = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        distance[frame] = bond_length(geometries[frame, atoms])[0]
    return distance


def get_angles(geometries, atoms):
    '''
    This functions wants an array with the geometries, and a list with a triplet of atoms
    '''
    number_of_steps = geometries.shape[0]
    number_of_atoms = geometries.shape[1]
    angle = np.empty(number_of_steps)
    for frame in range(number_of_steps):
        angle[frame] = bend_angle(geometries[frame, atoms])[0]
    return angle


def fit_gaussian(data, p0=[2, 0.1]):
    '''
    takes histogram and bins of same shape
    '''
    bins, hist = data[0], data[1]
    coeff, var_matrix = curve_fit(gaussian_distribution, bins, hist, p0=p0)
    # Get the fitted curve
    gauss_fit = gaussian_distribution(bins, *coeff)
    # get k of the oscillator that generates this distribution and append it to the coefficient list
    # sqrt(k/2pi kb T)
    temp = 330
    kb = boltzmann  # kb in hartree
    k = (temp * kb) / (coeff[1]**2)
    all_coefficients = np.append(coeff, k)
    # check if k replicates the real distribution
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
    Define the probability distribution function given K and T for testing purposes
    '''
    k, x0, temp = p
    kb = boltzmann
    return(1 / (np.sqrt(2 * np.pi * kb * temp / k)) * \
        np.exp(-(k / (2 * kb * temp)) * (x - x0)**2))


if __name__ == "__main__":
    msg = " angle_bond -p <path/to/trajectory>"
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
    args = parser.parse_args()
    main(args.p, args.st, args.et)
