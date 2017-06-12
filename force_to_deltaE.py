#!/usr/bin/env python

import sys
import numpy as np
from optparse import OptionParser

fname = sys.argv[1]

with open(fname) as f:
    content = f.readlines()

aantal_shots = len(content)/3

restraints=[]
forces=[]
stds=[]

for i in xrange(aantal_shots):
	lijn1=content[i*3].split()
	lijn2=content[i*3+1].split()
	lijn3=content[i*3+2].split()
	
	restraints.append(float(lijn1[2][:-1]))
	forces.append(float(lijn2[4]))
	stds.append(float(lijn3[4][:-1])**2)

energy = 0.0
std = 0.0
energy_list = []
energy_list.append(0.0)

for i in range(0,len(restraints)-1):
	energy+=((forces[i]+forces[i+1])/2.0)*(restraints[i+1]-restraints[i])
	energy_list.append(energy)
	std-=((stds[i]+stds[i+1])/2.0)*(restraints[i+1]-restraints[i])


std = np.sqrt(std)

print('Restraints:')
print restraints
print(' ')
print('Delta F = '+str(energy))
print('STD = '+str(std))
print('pKa = '+str(energy/(8.3145/1000*273.15*np.log(10))))
print('STD = '+str(std/(8.3145/1000*273.15*np.log(10))))
print(energy_list)

boltzmann = []

for i in range(len(energy_list)):
	boltzmann.append(np.exp(-energy_list[i]/(8.3145/1000*273.15)))

free_energy_prod = 0.0
free_energy_react = 0.0

cutoff = -1.5

for i in range(0,len(restraints)-1):
	if restraints[i] < cutoff:
		free_energy_prod += ((boltzmann[i]+boltzmann[i+1])/2.0)*(restraints[i+1]-restraints[i])
	else:
		free_energy_react += ((boltzmann[i]+boltzmann[i+1])/2.0)*(restraints[i+1]-restraints[i])

print('Delta F (Boltzmann): '+str(8.3145/1000*273.15*np.log(free_energy_react/free_energy_prod)))
print('pKa (Boltmann) = '+str(np.log(free_energy_react/free_energy_prod)/np.log(10)))
#	std_all=np.zeros(len(i))
# 	for j in xrange(len(i)):
# 		lijn=i[j].split()
# 		std_all[j] = lijn[1]
# 	std.append(np.max(std_all))
# 
# integral = 4.0/9.0*(deltaE_05)+5.0/18.0*(deltaE_01+deltaE_09)
# error =  4.0/9.0*(std[2])+5.0/18.0*(std[1]+std[2])

# print(' ')
# rint('Amount of MD-steps (min): '+str(np.min([len(contentf),len(contentg),len(contenth)])))
# print('Delta F = '+str(integral)+' eV.')
# print('Standard deviation = '+str(error)+' eV.')

