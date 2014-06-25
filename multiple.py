#!/usr/bin/env python

#
# multiple.py
#

# Import matplotlib and numpy modules
import matplotlib.pyplot as plt
import numpy as np

fig=plt.figure(1)

# Create the first subplot (first of 3 rows)
fig1=plt.subplot(311)

# Axes labels (vertical labels on the left)
plt.xlabel(  '$t$')
plt.ylabel('$S(t)$', rotation='horizontal')

# Load the data as arrays. We take the transpose because each array must be homogeneous
# The first (index 0) is always the time
redent = (np.loadtxt('RED-ENTROPY.dat')).T
cpent  = (np.loadtxt( 'CP-ENTROPY.dat')).T

# Plot the the graphs
rfig = plt.plot(redent[0], redent[1], color='red', 
        label='Redfield dynamics entropy')
cfig = plt.plot(cpent[0],  cpent[1],  color='blue',
        label='Completely positive dynamics entropy')

# Draw the grid
fig1.grid(True)

# Write the legend
plt.legend(('Redfield dynamics entropy',
    'Completely positive dynamics entropy'), loc='upper left',
    bbox_to_anchor=(0.2, 0.95))

# Create another label for the right vertical axis, leaving unaltered the horizontal one
ax = plt.twinx()
plt.ylabel('$\sigma(t)$', rotation='horizontal')

# Load the entropy production data into arrays
redprod = (np.loadtxt('RED-ENTROPY-PROD.dat')).T
cpprod  = (np.loadtxt('CP-ENTROPY-PROD.dat')).T

# Plot the graphs, draw the grid and write the legend
rpfig = plt.plot(redprod[0], redprod[1], 'y-', label='Redfield entropy production')
cpfig = plt.plot(cpprod[0], cpprod[1], 'c-', label='Completely positive entropy prod.')

plt.grid(True)
plt.legend(('Redfield entropy prod.', 'Completely posit. entropy prod.'),
        loc='upper left', bbox_to_anchor=(0.2, 0.6))




# Second subplot: time-evolution of the states in the Redfield dynamics
fig2 = plt.subplot(312)
redstat = (np.loadtxt('RED-EVOLUTION.dat')).T

# Every component of the Bloch vector must be plotted separetely
plt.plot(redstat[0],redstat[1])
plt.plot(redstat[0],redstat[2])
plt.plot(redstat[0],redstat[3])

fig2.grid(True)


# Third subplot: time-evolution of the states in the CP dynamics
fig3=plt.subplot(313)
cpstat = (np.loadtxt('CP-EVOLUTION.dat')).T

plt.plot(cpstat[0],cpstat[1])
plt.plot(cpstat[0],cpstat[2])
plt.plot(cpstat[0],cpstat[3])

fig3.grid(True)


plt.show()

