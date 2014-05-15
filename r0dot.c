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
 *       Filename:  r0dot.c
 *
 *    Description:  Check the sign of the time derivative at t=0
 *    			of the Bloch vector for an initial state with r=1 :
 *    			if negative, the dynamics is non positive.
 *
 *        Version:  1.0
 *        Created:  14/05/2014 18:07:11
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>
#include "funcs.h"


/* 
 *      FUNCTION  
 *         Name:  r0_dot
 *  Description:  Calculate the time derivative of r0 at time t=0 
 * 
 */
double r0_dot ( const gsl_matrix* K, const gsl_vector* R )
{
	/* time derivative of r at time t = 0 */
	double time_d ; 

	/* Take the spatial part (index = 1,2,3) of the Bloch generator */
	gsl_matrix_const_view L = gsl_matrix_const_submatrix ( K, 1, 1, 3, 3) ;

	/* spatial subcolumn of col. 0 */
	gsl_vector_const_view l = gsl_matrix_const_subcolumn ( K, 0, 1, 3 ) ;

	/* spatial part of Bloch vector R */
	gsl_vector_const_view r = gsl_vector_const_subvector ( R, 1, 3 ) ;

	/* L.r product */
	gsl_vector* y = gsl_vector_calloc(3) ;
	if ( gsl_blas_dgemv(CblasNoTrans, 1, &L.matrix, &r.vector, 0, y) != 0 )
		exit(EXIT_FAILURE) ;

	/* r.L.r scalar product */
	double scal1 ;
	if ( gsl_blas_ddot( &r.vector, y, &scal1 ) != 0 )
		exit(EXIT_FAILURE) ;

	/*  r.l scalar product */
	double scal2 ;
	if ( gsl_blas_ddot( &r.vector, &l.vector, &scal2 ) != 0 )
		exit(EXIT_FAILURE) ;

	time_d =  scal1 + scal2  ;
	 
	return (time_d) ;
}		/* -----  end of function r0_dot  ----- */
