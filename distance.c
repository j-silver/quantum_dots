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
 *       Filename:  distance.c
 *
 *    Description:  Distance in trace norm among the stationary and thermal
 *    			states.
 *
 *        Version:  1.0
 *        Created:  22/06/2014 22:21:50
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_vector.h>
#include "funcs.h"


/* 
 *      FUNCTION  
 *         Name:  dist
 *  Description:  distance between two states in Bloch vector form. It does not
 *  			check the positivity of the states
 * 
 */
double dist ( const gsl_vector* rho1, const gsl_vector* rho2 )
{
	double square = POW_2(VECTOR(rho1,1)-VECTOR(rho2,1)) +
			POW_2(VECTOR(rho1,2)-VECTOR(rho2,2)) +
			POW_2(VECTOR(rho1,3)-VECTOR(rho2.3)) ;
	double distance = gsl_sqrt(square) ;

	return (distance);
}		/* -----  end of function dist  ----- */


/* 
 *      FUNCTION  
 *         Name:  therm_state
 *  Description:  Return the thermal Gibbs state, given the parameter
 *  			beta and the Hamiltonian H
 * 
 */
gsl_vector therm_state ( const double beta, const gsl_matrix* H )
{
	/* Copy the Hamiltonian into a matrix object that can be destroyed */
	gsl_matrix* h = gsl_matrix_calloc(2, 2) ;
	gsl_matrix_memcpy ( h, H ) ;

	/* Find the eigenvalues.
	 * When written in Bloch form, states have an immediate spectral 
	 * decomposition:
	 *
	 * rho = (1+r)/2 (s0+r.s)/2 + (1-r)/2 (s0-r.s)/2
	 *
	 */
	double rsq = POW_2(VECTOR(
	return <+return_value+>;
}		/* -----  end of function therm_state  ----- */
