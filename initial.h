/* Copyright (c) 2014, Giuseppe Argentieri <giuseppe.argentieri@ts.infn.it>

 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*
 * 
 *
 *       Filename:  initial.h
 *
 *    Description:  Initial values
 *
 *        Version:  1.0
 *        Created:  16/05/2014 00:44:34
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_const_mksa.h>

#define	BOLTZ	GSL_CONST_MKSA_BOLTZMANN	/* Boltzmann constant  */
#define HBAR	GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR /* hbar */

const double omega_c = 1000 ;			/* critical ohmic frequency */
const double alpha = 5e-3 ;			/* coupling strength */
const double Delta = 8.0e+9 ;			/* pumping amplitude (GHz) */

double T = 0.1 ;				/* 
						 *  The real temperature is given by
						 *  temp = T*HBAR*Delta/BOLTZ
						 */
				
const double D = 1 ;                            /* normalized delta */
double OMEGA = 2 ;				/* normalized pumping frequency */

const double gamma0 = 0.05 ;			/* energy hopping between sites */
	
const double t_end = 200 ;			/* time end */
const double STEP = .01 ;			/* time step */

const double R[] = { 1, 0, -0.894, -0.447 } ;   /* initial state: |z,-> */
 

/* const double R[] = { 1, 0, 0.5, -0.4 } ; */ /* initial state with neg. e.p. */
/* const double r[] = { 1, 0, 1, 0 } ;  	 initial state with pos. t.d. */

