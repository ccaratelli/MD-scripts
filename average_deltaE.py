#!/usr/bin/env python

import sys
import numpy as np
from optparse import OptionParser

fname = sys.argv[1]
gname = sys.argv[2]

parser = OptionParser(
    """This script calculates the average vertical energy difference.""")
parser.add_option(
    "--start", "-s", default=0, type='int',
    help="Give the stepnumber from which to start the calculation.",
)

timestep = 0.5

(options, args) = parser.parse_args()

start = int(options.start)


def calc_cor(signal):
    def block_signal(signal):
        n = len(signal)
        if n % 2 == 1:
            signal = signal[:-1]
        return np.mean(signal.reshape(-1, 2), axis=1)

    n_entries = int(np.floor(np.log2(len(signal))))
    sigma_m = np.zeros(n_entries)
    sigma_m_std = np.zeros(n_entries)
#    for i in xrange(n_entries):
    for i in range(n_entries):
        n = float(len(signal))
        c0 = np.var(signal)
        sigma_m[i] = np.sqrt(c0 / (n - 1))
        sigma_m_std[i] = sigma_m[i] / np.sqrt(2 * (n - 1))
        signal = block_signal(signal)
    return sigma_m, sigma_m_std


with open(fname) as f:
    content = f.readlines()[start:]

g = open(gname + '.dat', 'w')
g2 = open(gname + '.block.dat', 'w')

g.write(
    '#Step' +
    '\t'
    'E1' +
    '\t\t' +
    'E2' +
    '\t\t'
    'Delta E (eV)' +
    '\t\t' +
    'Running avg' +
    '\n')

lijst = np.zeros(len(content))

steps = 10

counter = 0
delta = 0.0
delta_alles = 0.0

for i in range(len(content)):
    counter += 1
    line = content[i].split()
    E1 = float(line[3])
    E2 = float(line[4])
    delta = (E2 - E1) * 27.211385
    delta_alles += delta
    lijst[i] = delta
    g.write(str(i) + '\t' + str(E1) + '\t' + str(E2) + '\t' +
            str(delta) + '\t' + str(delta_alles / counter) + '\n')

sigma_m, sigma_m_std = calc_cor(lijst)

tijd = timestep * (2**np.arange(len(sigma_m)) - 1)

for j in range(len(sigma_m)):
    g2.write(str(tijd[j]) + '\t' + str(sigma_m[j]) +
             '\t' + str(sigma_m_std[j]) + '\n')

print('Total average delta E = ' + str(delta_alles / counter) + ' eV.')
print('Sample standard deviation = ' + str(max(sigma_m)) + ' eV.')

f.close
g.close
g2.close
