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
 *       Filename:  entropy.c
 *
 *    Description:  Calculate the entropy production for a given state, generator
 *    			and equilibrium state
 *
 *        Version:  1.0
 *        Created:  05/05/2014 20:53:01
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "funcs.h"

/* 
 *      FUNCTION  
 *         Name:  entropy_production
 *  Description:  
 * 
 */
double entropy_production ( const gsl_vector* rho, const gsl_vector* rhoeq,
		const gsl_matrix* L  )
{
	/* l1, l2, l3 */
	double l[3] ; unsigned int i, j ;
	for ( i = 1 ; i < 3 ; i++ )
	{
		l[i] = 0 ;
		for ( j = 0 ; j < 3 ; j++ )
			l[i] += gsl_matrix_get(L,i,j)*VECTOR(rho,j) ;
	}	

	/* L[rho] */
	double Lr = 0 ;
	for ( i = 1 ; i < 3 ; i++ )
		Lr += l[i]*VECTOR(rho,i) ;

	/* L[rhoeq] */
	double Leq = 0 ;
	for ( i = 1 ; i < 3 ; i++ )
		Leq += l[i]*VECTOR(rhoeq,i) ;

	/* r , req */
	double r, req ;
	r = req = 0 ;

	r = gsl_hypot3(VECTOR(rho,1),VECTOR(rho,2),VECTOR(rho,3)) ;
	req = gsl_hypot3(VECTOR(rhoeq,1),VECTOR(rhoeq,2),
		VECTOR(rhoeq,3)) ;

	/* internal entropy s */
	double s ;
	if ( r < 1 && req < 1 )
		s = -(gsl_sf_log((1+r)/(1-r))*Lr/r - gsl_sf_log((1+req)/(1-req))*Leq/req) ;
	else 
		s = 0 ;

	return s;
}		/* -----  end of function entropy_production  ----- */


/* 
 *      FUNCTION  
 *         Name:  entropy_of_state
 *  Description:  Calculate the Von Neumann entropy of state 'rho'
 * 
 */
double entropy_of_state ( const gsl_vector* rho )
{
	double entr = 0 ;

	/* Finding the eigenvalues */
	gsl_eigen_herm_workspace* rho_ei = gsl_eigen_herm_alloc(2);
	gsl_matrix_complex* dens = gsl_matrix_complex_calloc (2,2);
	gsl_matrix_complex_set (dens, 0, 0, gsl_complex_rect(1+VECTOR(rho, 3),0));
	gsl_matrix_complex_set (dens,0,1,gsl_complex_rect(VECTOR(rho,1),-VECTOR(rho,2)));
	gsl_matrix_complex_set (dens,1,0,gsl_complex_rect(VECTOR(rho,1),VECTOR(rho,2)));
	gsl_matrix_complex_set (dens,1,1,gsl_complex_rect(1-VECTOR(rho,3),0));
	gsl_matrix_complex_scale (dens, gsl_complex_rect(0.5,0));
	gsl_vector* eigenvalues = gsl_vector_calloc(2) ;
	gsl_eigen_herm (dens, eigenvalues, rho_ei) ;
	
	/* Calculating entropy */
	entr = - (VECTOR(eigenvalues,0)*gsl_sf_log(VECTOR(eigenvalues,0)) +
			VECTOR(eigenvalues,1)*gsl_sf_log(VECTOR(eigenvalues,1))) ;

	return (entr);
}		/* -----  end of function entropy_of_state  ----- */
