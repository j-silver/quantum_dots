#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()

plt.grid(True)

plt.xlabel('$t \Delta$', fontsize='large')
plt.ylabel('$r$', rotation='horizontal', fontsize='large')

STATE = (np.loadtxt('RED-EVOLUTION_200.dat')).T

plt.plot(STATE[0], STATE[4])

plt.draw()
plt.show()
