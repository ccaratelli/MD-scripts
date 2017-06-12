#!/usr/bin/env python

import sys
import numpy as np
from optparse import OptionParser

with open("lambda0.1.dat") as f:
    contentf = f.readlines()
with open("lambda0.5.dat") as g:
    contentg = g.readlines()
with open("lambda0.9.dat") as h:
    contenth = h.readlines()

with open("lambda0.1.block.dat") as f2:
    contentf2 = f2.readlines()
with open("lambda0.5.block.dat") as g2:
    contentg2 = g2.readlines()
with open("lambda0.9.block.dat") as h2:
    contenth2 = h2.readlines()

lijn = contentf[-1].split()
deltaE_01 = float(lijn[-1])

lijn = contentg[-1].split()
deltaE_05 = float(lijn[-1])

lijn = contenth[-1].split()
deltaE_09 = float(lijn[-1])

std=[]

for i in [contentf2,contentg2,contenth2]:
	std_all=np.zeros(len(i))
	for j in xrange(len(i)):
		lijn=i[j].split()
		std_all[j] = lijn[1]
	std.append(np.max(std_all))

integral = 4.0/9.0*(deltaE_05)+5.0/18.0*(deltaE_01+deltaE_09)
error =  4.0/9.0*(std[2])+5.0/18.0*(std[1]+std[2])

print(' ')
print('Amount of MD-steps (min): '+str(np.min([len(contentf),len(contentg),len(contenth)])))
print('Delta F = '+str(integral)+' eV.')
print('Standard deviation = '+str(error)+' eV.')

f.close
g.close
h.close

f2.close
g2.close
h2.close

