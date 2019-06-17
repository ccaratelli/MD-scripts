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
    Loads molecular geometries, generates file with time evolution
    of given angles and bonds. Prints force parameters for a 
    harmonic oscillator that reproduces the behavior of 
    each bond/angle to use as input for biased simulations
    """

    # Timestep in fs
    timestep = 0.5

    # Store the geometries from the trajectory as a numpy array
    xyz_file = XYZFile(file_name)
    geometries = xyz_file.geometries[start_step:end_step]

    # Read atom list from input file (input is 1-based indexing)
    with open('atoms.txt', 'r') as f:
        atoms_input = json.load(f)
    atoms = [np.array(a) - 1 for a in atoms_input]

    # Calculate bonds and angles
    time = (np.arange(geometries.shape[0]) + 1) * timestep
    bonds_angles = [get_bonds_angles(geometries, i) for i in atoms]
    labels = [convert_label(i) for i in atoms_input]
    units = ['Angle (rad)' if len(
        i) == 3 else 'Bond (a.u.)' for i in atoms_input]

    # Compute histograms and saves results
    for i, qty in enumerate(bonds_angles):
        all_distr, coefficients = generate_histogram(qty, temp)
        name = labels[i]
        unit = units[i]
        np.savetxt(name + "-hist.dat", all_distr)
        np.savetxt(name + "-time.dat", np.stack((time, qty)).transpose())
        np.savetxt(name + "-coeff.dat", coefficients, fmt='%1.3f')
        plot_all(all_distr, qty, coefficients, name, time, unit)

    # Store in a pandas dataframe for further analysis
    all_data = pd.DataFrame(
        data=np.stack(bonds_angles).transpose(),
        index=time,
        columns=labels)
    all_data.to_csv("all_data.dat", sep='\t')


def generate_histogram(colvar, temp):
    """
    Calculates a histogram from a quantity evolution during the simulation
    """
    # Define the histogram and shift the bins to have data on the centres
    hist, bin_edges = np.histogram(colvar, bins=50, density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    histogram = np.stack((bin_centres, hist))
    gaussian_distr, oscillator_distr, coefficients = fit_distribution(histogram, temp)
    all_distr = np.stack((bin_centres, hist, gaussian_distr, oscillator_distr), axis=1)
    return all_distr, coefficients


def fit_distribution(data, temp, p0=[2, 0.1]):
    """
    Takes histogram and bins of same shape and returns a fitted gaussian distribution
    and the force constants of the corresponding harmonic oscillator
    """
    bins, hist = data[0], data[1]
    coeff, var_matrix = curve_fit(gaussian_distribution, bins, hist, p0=p0)

    # Fit a gaussian distribution to the bond/angle distribution
    # We do not fit directly the oscillator distribution because sometimes the fit fails
    gauss_fit = gaussian_distribution(bins, *coeff)

    # Obtain force constant for the oscillator that generates this distribution
    # append it to the coefficient list.
    # sigma = sqrt(k/2pi kb T)
    kb = boltzmann  # kb in hartree
    k = (temp * kb) / (coeff[1]**2)
    all_coefficients = np.append(coeff, k)

    # Check if k replicates the real distribution
    coeff_distr = k, coeff[0], temp
    osc_fit = oscillator_distribution(bins, *coeff_distr)
    return gauss_fit, osc_fit, all_coefficients


def gaussian_distribution(x, *p):
    """
    Defines a gaussian function to fit data. A = 1 / np.sqrt(2 * np.pi * sigma**2)
    """
    mu, sigma = p
    return (1 / (np.sqrt(2 * np.pi) * sigma)) * \
        np.exp(-(x - mu)**2 / (2. * sigma**2))


def oscillator_distribution(x, *p):
    """
    Define the probability distribution function at given k and T for testing purposes.
    This distribution should be the same as the gaussian one.
    """
    k, x0, temp = p
    kb = boltzmann
    return(1 / (np.sqrt(2 * np.pi * kb * temp / k)) *
           np.exp(-(k / (2 * kb * temp)) * (x - x0)**2))


def get_bonds_angles(geometries, atoms):
    """
    This functions takes an array with the geometries, and a list with a group of atoms
    and returns the bond or angle evolution during the simulation

    """
    number_of_steps, number_of_atoms, _ = geometries.shape
    colvar = np.empty(number_of_steps)
    if len(atoms) == 2:
        for frame in range(number_of_steps):
            colvar[frame] = bond_length(geometries[frame, atoms])[0]
    elif len(atoms) == 3:
        for frame in range(number_of_steps):
            colvar[frame] = bend_angle(geometries[frame, atoms])[0]
    return colvar


def convert_label(colvar):
    """
    Generates a label with the list of atoms  [70, 170] --> Bond_70_170
    """
    if len(colvar) == 3:
        label = "angle_"
    else:
        label = "bond_"
    lab = '_'.join(str(x) for x in colvar)
    return label + lab


def plot_all(all_distr, qty, coefficients, name, time, unit):
    """
    Plots all data
    """
    fig, (p1, p2) = plt.subplots(1, 2, figsize=(
        12, 3), gridspec_kw={'width_ratios': [3, 1]})

    # Plot with the time evolution of the bond/length
    p1.set_xlabel('Time (ps)')
    p1.set_ylabel(unit)
    p1.set_title('Evolution of ' + name)
    p1.plot(time * 0.001, qty)

    # Plot with the distribution + fit
    p2.set_xlabel('Distribution')
    p2.axes.get_yaxis().set_visible(False)
    p2.plot(all_distr[:, 1], all_distr[:, 0])
    p2.plot(all_distr[:, 2], all_distr[:, 0])

    # Annotate the values for the distribution
    textstr = '\n'.join((
        r'$\mu=%.2f$' % (coefficients[0]),
        r'$\sigma=%.2f$' % (coefficients[1]),
        r'$k=%.2f$' % (coefficients[2])))
    p2.text(0.7, 0.95, textstr, transform=p2.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(facecolor='orange', alpha=0.7))

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)
    plt.savefig(name + ".png")


if __name__ == "__main__":
    msg = "angle_bond -p <path/to/trajectory> -st <start frame> -et <end frame>  -t <temperature>"
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

