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
 *       Filename:  station.c
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
int stationary ( const gsl_matrix* M, gsl_vector* X )
{
	/* Copy the dissipator matrix in a temporary local matrix m */
	gsl_matrix* m = gsl_matrix_calloc ( 4, 4 ) ;
	gsl_matrix_memcpy ( m, M ) ;
	gsl_vector_view x = gsl_vector_subvector ( X, 1, 3 ) ;

	/* Create a submatrix view of the spatial part of m and a vector view
	 * of the spatial part of the 0-th column, which goes into -b in the system
	 * A X = b */
	gsl_matrix_view A = gsl_matrix_submatrix ( m, 1, 1, 3, 3 ) ;
	gsl_vector_view b = gsl_matrix_subcolumn ( m, 0, 1, 3 ) ;
	int status1 = gsl_vector_scale ( &b.vector, -1.0 ) ;	

	/* Solve the system A x = b using Householder transformations.
	 * Changing the view x of X => also X is changed, in the spatial part */
	int status2 = gsl_linalg_HH_solve ( &A.matrix, &b.vector, &x.vector ) ;

	return status1 + status2 ;
}		/* -----  end of function stationary  ----- */
