#!/usr/bin/env python

#
# entp_plots.py
#

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

# Reading data
ent = (np.loadtxt('t-fixed-2.0.DAT')).T

plt.xlabel('$t \Delta$')
plt.ylabel('$\sigma(t)$', rotation='horizontal')

plt.plot(ent[0], ent[1], 'y-')

plt.grid(True)

plt.show()
