#!/usr/bin/env python

#
# Positivity.py
#

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm

# Load data
PTS = (np.loadtxt ('POS_VIOLATIONS')).T

X = PTS[0]
Y = PTS[1]
Z = PTS[2]

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

surf = ax.scatter( X, Y, Z, s=0.1, color='r')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
