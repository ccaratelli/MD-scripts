#!/usr/bin/env python

import sys
import numpy as np
from optparse import OptionParser

fname = sys.argv[1]
gname = sys.argv[2]

parser = OptionParser("""This script calculates the average force.""")
parser.add_option(
    "--start", "-s", default=0, type='int',
    help="Give the stepnumber from which to start the calculation.",
)

(options, args) = parser.parse_args()

start=int(options.start)

def calc_cor(signal):
    """ This function gives sigma from a normal distribution """
    def block_signal(signal):
        n = len(signal)
        if n%2 == 1:
            signal = signal[:-1]
        return np.mean(signal.reshape(-1,2), axis=1)

    n_entries = int(np.floor(np.log2(len(signal))))
    sigma_m = np.zeros(n_entries)
    sigma_m_std = np.zeros(n_entries)
    for i in xrange(n_entries):
        n = float(len(signal))
        c0 = np.var(signal)
        sigma_m[i] = np.sqrt(c0/(n-1))
        sigma_m_std[i] = sigma_m[i]/np.sqrt(2*(n-1))
        signal = block_signal(signal)
    return sigma_m, sigma_m_std

with open(fname) as f:
    content = f.readlines()[start:]

g = open(gname+'.dat', 'w')
g2 = open(gname+'.block.dat', 'w')

g.write('#Step'+'\t''Force'+'\t\t'+'Running avg'+'\n')

counter=0
kracht = 0.0
kracht_alles = 0.0

constant = 1500.0
set_value = 0.351880

lijst = np.zeros(len(content)-1)

for i in range(1,len(content)):
	counter += 1
	line = content[i].split()
	kracht = -constant*(float(line[3])-set_value)
	kracht_alles += kracht
	lijst[i-1] = kracht
	g.write(str(i)+'\t'+str(kracht)+'\t'+str(kracht_alles/counter)+'\n')

sigma_m, sigma_m_std = calc_cor(lijst)

tijd = 2**np.arange(len(sigma_m))-1

for j in xrange(len(sigma_m)):
        g2.write(str(tijd[j])+'\t'+str(sigma_m[j])+'\t'+str(sigma_m_std[j])+'\n')

print('Restraint at '+str(set_value)+'.')
print('Total average Force = '+str(str(kracht_alles/counter))+' kJ/mol.')
print('Sample standard deviation = '+str(max(sigma_m))+' kJ/mol.')

f.close
g.close
g2.close

