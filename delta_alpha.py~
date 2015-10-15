#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

#plt.grid(True)

plt.xlabel('$\\alpha$', fontsize='large')
plt.ylabel('$\\frac{\Delta(\\alpha)}{(1+\\alpha^2)^2}$', rotation='horizontal', fontsize='large')

DATA = (np.loadtxt('D_ALPHA.dat')).T

plt.plot(DATA[0], DATA[1])

#plt.text(6, -0.02, '$\\frac{\kappa_B T}{\hbar\Delta}=1$\n$\\lambda=0.005$', bbox={'facecolor':'white','alpha':0.8})

plt.draw()
plt.show()
