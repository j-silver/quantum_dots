# Thermodynamics of open quantum systems

## Physical context

This repository contains the implementation of a numerical model of an open quantum system studied in [1] and [2].

The articles analyse a simple micro-circuit made of a 3-site quantum dot in which three electrons are driven by a periodic potential while the circuit is immersed in a noisy thermal bath.

The modelisation leads to a dynamics governed by a typical Lindbladian Master Equation which accounts for both the driving Hamiltonian and the environment dissipative effect.

Depending on the adopted approximation technique the resulting dynamics can be either *completely positive* (CP) or non-CP. 

I developed my simulation in order to track the time-evolution of the *internal entropy production* and of the *global state* of the system. The former allows to study possible violation of the 2^n^d law of thermodynamics (highlighted by negative values of the entropy production), while the time-evolution of the state contains all the relevant information of the system. In particular it is possible to determine the asymptotic current at regime and its dependence on the physical parameters at play (driving frequency, amplitude and so on).


## Code implementation

The model is implemented in C and the code is split into many files, each of them devoted to the calculation of a specific quantity. 

In particular:

* the files Img... and Reg... are used to evaluate the imaginary and real part of the integrals (85)-(98) appearing in [1], which are used to build the Master Equations. 

* the files red\_gen.c, cp\_gen.c implement the matrix generators of the Redfield (non-CP) and CP Master Equations.

* the files red\_evol.c, cp\_evol.c solve the differential Master Equation in the both cases, by adopting standard Runge-Kutta methods. They write the time-evolution of the *entropy production* and of the *Bloch vector* into two files, which are then used to plot the results by usnig the various Python scripts.

* the other C files are used to evaluate other quantities such as the asymptotic current (asymptotic.c), the current's' dependence on some parameters like the frequency Omega or the temperature (current\_omegad.c, current\_t.c), for a 3D visualization of the Bloch sphere where the initial entropy violations occur (3Dposit.c) or for test purposes.


The integrals, the matrix manipulation and the differential equations are computed by using specific functions from the __GNU Scientific Library GSL__ [link]https://www.gnu.org/software/gsl/[link].






 




[1] [link]http://link.springer.com/article/10.1007/s10955-015-1210-4 
	(arXiv version: [link]https://arxiv.org/abs/1502.00864)

[2] [link]http://iopscience.iop.org/article/10.1209/0295-5075/107/50007/meta
	(arXiv version: [link]https://arxiv.org/abs/1408.4589)


