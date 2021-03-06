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
 *       Filename:  stationary.c
 *
 *    Description:  Find the stationary state for the given generator
 *
 *        Version:  1.0
 *        Created:  05/05/2014 15:27:48
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "funcs.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


/* 
 *      FUNCTION  
 *         Name:  stationary
 *  Description:  Given the dissipator in Bloch form, reduce to a 3x3 problem and store
 *  			the stationary state in the 3x1 vector *X 
 * 			
 * 			M X = 0
 * 
 * 	        	|  0    0    0    0  |  | 1  |    0
 * 		        | M10  M11  M12  M13 |  | X1 |    0
 * 		        | M20  M21  M22  M23 |  | X2 | =  0
 * 			| M30  M31  M32  M33 |  | X3 |    0
 *
 *
 * 			A x = b
 *
 * 			| M11  M12  M13 |  | X1 |   | -M10 |
 * 			| M21  M22  M23 |  | X2 | = | -M20 |
 * 			| M31  M32  M33 |  | X3 |   | -M30 |
 */
int stationary ( const gsl_matrix* M, gsl_vector* stat_state )
{
	/* Store space for the stationary state */
	gsl_vector* req = gsl_vector_calloc ( 4 ) ;
	gsl_vector_set ( req, 0, 1 ) ;

	/* Copy the dissipator matrix in a temporary local matrix m
	 * (because the algorithm destroys it...) */
	gsl_matrix* m = gsl_matrix_calloc ( 4, 4 ) ;
	gsl_matrix_memcpy ( m, M ) ;

	/* Create a view of the spatial part of vector req */
	gsl_vector_view x = gsl_vector_subvector ( req, 1, 3 ) ;

	/* Create a submatrix view of the spatial part of m and a vector view
	 * of the spatial part of the 0-th column, which goes into -b in the system
	 * A x = b */
	gsl_matrix_view A = gsl_matrix_submatrix ( m, 1, 1, 3, 3 ) ;
	gsl_vector_view b = gsl_matrix_subcolumn ( m, 0, 1, 3 ) ;
	int status1 = gsl_vector_scale ( &b.vector, -1.0 ) ;	

	/* Solve the system A x = b using Householder transformations.
	 * Changing the view x of req => also req is changed, in the spatial part */
	int status2 = gsl_linalg_HH_solve ( &A.matrix, &b.vector, &x.vector ) ;

	/* Set the returning value for the state stat_state */
	*stat_state = *req ;

	/* Free memory */
	gsl_matrix_free(m) ;
	
	return status1 + status2 ;
}		/* -----  end of function stationary  ----- */



/* 
 *      FUNCTION  
 *         Name:  J
 *  Description:  Ohmic spectral density (divided by alpha)
 * 
 */
double J ( double w, double oc )
{
	double W = w*exp(-w/oc) ;
	return (W);
}		/* -----  end of function J  ----- */


/* 
 *      FUNCTION  
 *         Name:  polarization
 *  Description:  Polarization P through the CP formula 
 * 
 */
double polarization ( void* params, double oc )
{
	struct f_params* pars = (struct f_params*) params ;
	double Omega, omega_1, b ;
	Omega   = pars->Omega ;
	omega_1 = pars->omega_1 ;
	b       = pars->beta ;

	double P, num, den, add1, add2, Jplus, Jminus ;
	Jplus  = J(omega_1 + Omega, oc) ;
	Jminus = J(omega_1 - Omega, oc) ;
	add1 = POW_2(omega_1 - Omega)*Jplus ;
	add2 = POW_2(omega_1 + Omega)*Jminus ;
	num = add1 + add2 ;
	den = add1/(tanh((omega_1+Omega)*b/2)) +
		add2/(tanh((omega_1-Omega)*b/2)) ;
	P = num/den ;

	return (P);
}		/* -----  end of function polarization  ----- */
