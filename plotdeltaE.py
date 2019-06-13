#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib  # comment for figure
matplotlib.use('TkAgg')  # comment for figure

# output="deltae.png" #decomment for figure
l1 = np.loadtxt("lambda0.1.dat").transpose()
l5 = np.loadtxt("lambda0.5.dat").transpose()
l9 = np.loadtxt("lambda0.9.dat").transpose()

l1[0] = l1[0] * 0.5
l5[0] = l5[0] * 0.5
l9[0] = l9[0] * 0.5
font_size = 18
plt.xlabel("Time (fs)", fontsize=font_size)
plt.ylabel(r"$\Delta$E (eV)", fontsize=font_size)

plt.plot(l1[0], l1[3], l5[0], l5[3], l9[0], l9[3], color='0.75')
plt.plot(l1[0], l1[4], label=r'$\eta$ = 0.1')
plt.plot(l5[0], l5[4], label=r'$\eta$ = 0.5')
plt.plot(l9[0], l9[4], label=r'$\eta$ = 0.9')

plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.show(block=True)
# plt.savefig(output, format='png', dpi=300)  #decomment for figure
