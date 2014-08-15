#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np


#ICP  = 0.894427191
#IRED = 0.877812926

#t = np.arange( 0.00, 600.00, .01 )

#RED = IRED*np.ones_like(t)
#CP  = ICP *np.ones_like(t)

#plt.plot(t, CP, 'g--')
#plt.plot(t, RED, 'r--')

RCUR = (np.loadtxt('RED-CURRENT.dat')).T
CCUR = (np.loadtxt('CP-CURRENT.dat')).T

plt.plot(RCUR[0], RCUR[1], 'r', label='Redfield dynamics current')
plt.plot(CCUR[0], CCUR[1], 'g', label='CP dynamics current')

plt.xlabel('$t \Delta$')
plt.ylabel('$I/I_0$', rotation='horizontal')

plt.text(150, 0.6, '$\Omega/\Delta=2$\n$\kappa_B T/\hbar \Delta=10$\n$\\alpha=0.005$\n$\\rho_0=\\vert z; -\\rangle$', bbox={'facecolor':'white'})
#\{1, 0, 0.5, -0.4\}$', bbox={'facecolor':'white'}) #
plt.legend(('Redfield dynamics current', 'CP dynamics current'), loc='lower right',
		bbox_to_anchor=(0.95, 0.5))
plt.grid(True)

plt.draw()
plt.show()

