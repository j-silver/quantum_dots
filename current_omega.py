#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
#fig.suptitle('Asymptotic current as a function of the frequency $\Omega$')

plt.xlabel('$\Omega/\\Delta$')
plt.ylabel('$I/I_0$')

plt.xscale('log')
plt.grid('on')

#plt.text(0.005, 0.15, '$\\frac{\kappa_B T}{\hbar\Delta}=1$\n$\\alpha=0.005$', bbox={'facecolor':'white','alpha':0.8})

R_CURRT = (np.loadtxt('RED-STAT-CURR-O.dat')).T
C_CURRT = (np.loadtxt('CP-STAT-CURR-O.dat')).T


rfig = plt.plot(R_CURRT[0], R_CURRT[1], color='red', label='Redfield dynamics')
cfig = plt.plot(C_CURRT[0], C_CURRT[1], color='blue', label='CP dynamics')

plt.legend(('Redfield dynamics', 'CP dynamics'))

plt.show()
