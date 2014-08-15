#!/usr/bin/env python

#
# entropy_plot.py
#

import matplotlib.pyplot as plt
import numpy as np
# import csv

fig = plt.figure()


# Reading data
redent = (np.loadtxt('RED-ENTROPY.dat')).T
cpent  = (np.loadtxt( 'CP-ENTROPY.dat')).T

redprod = (np.loadtxt('RED-ENTROPY-PROD.dat')).T
cpprod  = (np.loadtxt('CP-ENTROPY-PROD.dat')).T

#
# Subplot n.1 : Entropy evolution
#
plt.subplot(211)
plt.title('Entropy time-evolution')

# Setup labels
plt.xlabel('$t \Delta$')
plt.ylabel('$S(t)$', rotation='horizontal')

#
# Text box
#
plt.text(100, 0.1, '$\Omega/\Delta=2$\n$\kappa_B T/\hbar \Delta=0.1$\n$\\alpha=0.005$\n$\\rho_0=\\vert z;-\\rangle$', bbox={'facecolor':'white'}) 
#plt.text(250, 0.5, '$\Omega/\Delta=2$\n$\kappa_B T/\hbar \Delta=0.1$\n$\\alpha=0.005$\n$\\rho_0=\{1, 0, 0.5, -0.4\}$', bbox={'facecolor':'white'})


# Plotting
rfig = plt.plot(redent[0], redent[1], color='red', 
		label='Redfield dynamics entropy')
cfig = plt.plot(cpent[0],  cpent[1],  color='blue',
		label='Completely positive dynamics entropy')

# Maximum entropy 
maxent = np.log(2.0)*np.ones_like(redent[0])
plt.plot(redent[0], maxent)

plt.grid(True)

plt.legend(('Redfield dynamics entropy',
	'Completely positive dynamics entropy'), loc='upper left',
	bbox_to_anchor=(0.2, 0.95))

#ax = plt.twinx()

#
# Subplot n.2 : Entropy production
#
plt.subplot(212)
plt.title('Internal entropy production')

plt.xlabel('$t \Delta$')
plt.ylabel('$\sigma(t)$', rotation='horizontal')

rpfig = plt.plot(redprod[0], redprod[1], 'y-', label='Redfield entropy prod.')
cpfig = plt.plot(cpprod[0], cpprod[1], 'c-', label='Completely posit. entropy prod.')

plt.grid(True)

plt.legend(('Redfield entropy prod.', 'Completely posit. entropy prod.'),
		loc='upper left', bbox_to_anchor=(0.2, 0.8))

plt.show()

