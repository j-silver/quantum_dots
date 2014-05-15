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
 *
 *
 *       Filename:  evol.c
 *
 *    Description:  Evolution of the density matrix with a given generator
 *
 *        Version:  1.0
 *        Created:  02/05/2014 14:40:25
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "funcs.h"


/* 
 *      FUNCTION  
 *         Name:  generator
 *  Description:  Setting the dydt in Bloch form
 * 
 */
int generator ( double t, const double y[], double dydt[], void* PARS )
{
	/* Extracting the matrix address from PARS */
	gsl_matrix* M = (gsl_matrix*) PARS ;

	/* Generator */
	unsigned int i , j ;
	for ( i = 0; i < 4 ; i++ )
		dydt[i] = 0 ;

	for ( i = 1 ; i < 4 ; i++ )
		for ( j = 0 ; j < 4 ; j++ )
			dydt[i] = dydt[i] + gsl_matrix_get( M, i, j )*y[j] ;

	return GSL_SUCCESS;
}		/* -----  end of function generator  ----- */


/* 
 *      FUNCTION  
 *         Name:  jac
 *  Description:  Jacobian matrix of the generator (which is exactly the Bloch matrix)
 * 
 */
int jac ( double t, const double y[], double dfdy[] , double dfdt[], void* PARS )
{
	/* Taking the generator Bloch's matrix address from PARS */
	gsl_matrix* bloch = (gsl_matrix*) PARS ;

	/*  Creating a matrix view of the array dfdy */
	gsl_matrix_view m = gsl_matrix_view_array ( dfdy , 4, 4 ) ; 	
	
	/* Initializing the jacobian matrix with the Bloch generator */
	unsigned int i, j ; 
	for ( i = 0 ; i < 4 ; i++ ) 
		for ( j = 0 ; j < 4 ; j++ ) 
			gsl_matrix_set ( &m.matrix, i, j, gsl_matrix_get(bloch,i,j) ) ;

	for ( i = 0 ; i < 4 ; i++ )
		dfdt[i] = 0 ;

	return GSL_SUCCESS;
}		/* -----  end of function jac  ----- */


/* 
 *      FUNCTION  
 *         Name:  evol
 *  Description:  
 * 
 */
int evol ( double t, gsl_vector* state, double step,
		gsl_odeiv2_evolve* e, gsl_odeiv2_control* c, gsl_odeiv2_step* s,
		gsl_odeiv2_system* sys )
{
	/* Creating the array 'rho' pointing to vector 'state' */
	double rho[4] ;
	unsigned int i ;
	for ( i = 0 ; i < 4 ; i++ ) 
		rho[i] = VECTOR( state, i ) ;

	/*  evolving the array 'rho'  */
	int status = gsl_odeiv2_evolve_apply_fixed_step	( e, c, s, sys, &t, step, rho ) ;
	if ( status != GSL_SUCCESS )
	{
		printf("STATUS: %d\n", status) ;
		exit(1) ;
	}

	/* updating vector 'state' */
	for ( i = 0 ; i < 4 ; i++ )
		gsl_vector_set( state, i, rho[i] ) ;

	return status ;
}		/* -----  end of function evol  ----- */



