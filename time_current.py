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

plt.text(200, 0.4, '$\Omega/\Delta=2$\n$T=0.1$\n$\\alpha=0.005$\n$\\rho_0=\{1, 0, 0.5, -0.4\}$', bbox={'facecolor':'white'}) #$\\vert z; -\\rangle$', bbox={'facecolor':'white'}

plt.legend(('Redfield dynamics current', 'CP dynamics current'), loc='lower right',
		bbox_to_anchor=(0.95, 0.05))
plt.grid(True)

plt.draw()
plt.show()

