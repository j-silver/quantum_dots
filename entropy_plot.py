#!/usr/bin/env python

#
# entropy_plot.py
#

import matplotlib.pyplot as plt
import numpy as np
import csv

fig = plt.figure()
fig.suptitle('Entropy time-evolution')

plt.xlabel(  '$t \Delta$')
plt.ylabel('$S(t)$', rotation='horizontal')


#plt.text(350, .6, first, bbox={'facecolor':'white','alpha':0.5})
plt.text(250, 0.5, '$\Omega/\Delta=2$\n$T=0.006$\n$\\alpha=0.005$\n$\\rho_0=\\vert z;-\\rangle$', bbox={'facecolor':'white'}) #\{1, 0, 0.5, -0.4\}$', bbox={'facecolor':'white'})

redent = (np.loadtxt('RED-ENTROPY.dat')).T
cpent  = (np.loadtxt( 'CP-ENTROPY.dat')).T

rfig = plt.plot(redent[0], redent[1], color='red', 
		label='Redfield dynamics entropy')
cfig = plt.plot(cpent[0],  cpent[1],  color='blue',
		label='Completely positive dynamics entropy')

maxent = np.log(2.0)*np.ones_like(redent[0])
plt.plot(redent[0], maxent)

plt.grid(True)

plt.legend(('Redfield dynamics entropy',
	'Completely positive dynamics entropy'), loc='upper left',
	bbox_to_anchor=(0.2, 0.95))


ax = plt.twinx()

plt.ylabel('$\sigma(t)$', rotation='horizontal')

redprod = (np.loadtxt('RED-ENTROPY-PROD.dat')).T
cpprod  = (np.loadtxt('CP-ENTROPY-PROD.dat')).T

rpfig = plt.plot(redprod[0], redprod[1], 'y-', label='Redfield entropy prod.')
cpfig = plt.plot(cpprod[0], cpprod[1], 'c-', label='Completely posit. entropy prod.')

plt.grid(True)


plt.legend(('Redfield entropy prod.', 'Completely posit. entropy prod.'),
		loc='upper left', bbox_to_anchor=(0.2, 0.8))

plt.show()

