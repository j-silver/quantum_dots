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
 *       Filename:  posit.c
 *
 *    Description:  Determine the positivity violation of the Redfield dynamics
 *    			on a given initial pure state
 *
 *        Version:  1.0
 *        Created:  20/07/2014 09:59:35
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "initial.h"
#include "funcs.h"

#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	double o_eff = gsl_hypot(OMEGA, D);

	/* sigma^z eigenstates, Bloch vectorial form */
	double z_minus[] = { 1, 0, -OMEGA/o_eff, -D/o_eff };
	double z_plus[]  = { 1, 0,  OMEGA/o_eff,  D/o_eff };
	gsl_vector* z_m = gsl_vector_calloc(4);
	gsl_vector* z_p = gsl_vector_calloc(4);

	int i;                                  
	for ( i = 0; i < 4; i++ ) {
		gsl_vector_set( z_m, i, z_minus[i] );
		gsl_vector_set( z_p, i, z_plus[i]  );
	}

	/* parameters */
	double r, theta;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */


/* 
 *      FUNCTION  
 *         Name:  state_matrix
 *  Description:  return the density patrix of the pure state |psi><psi|
 * 
 */
gsl_matrix_complex* state_matrix ( double r, double theta )
{
	gsl_matrix_complex* psi = gsl_matrix_complex_calloc (2, 2);
	gsl_matrix_complex_set ( psi, 1, 1, gsl_complex_rect(r*r, 0) );
	gsl_matrix_complex_set ( psi, 1, 2, gsl_complex_rect(r*sqrt(1-r*r)*cos(theta),
				r*sqrt(1-r*r)*sin(theta)) );
	gsl_matrix_complex_set ( psi, 2, 1, gsl_complex_rect(r*sqrt(1-r*r)*cos(theta),
				-r*sqrt(1-r*r)*sin(theta)) );
	gsl_matrix_complex_set ( psi, 2, 2, gsl_complex_rect(1-r*r, 0) );

	return psi;
}		/* -----  end of function state_matrix  ----- */


/* 
 *      FUNCTION  
 *         Name:  rotate
 *  Description:  return the density matrix in the rotated frame of reference
 *  			of an angle OMEGA*t
 * 
 */
gsl_matrix_complex* rotate ( gsl_matrix_complex* m, double t )
{
	/* copy m into a working matrix M */
	gsl_matrix_complex* M = gsl_matrix_complex_calloc ( 2, 2);
	gsl_matrix_complex_memcpy( M, m );

	/* rotation matrix */
	gsl_matrix_complex* Rot = gsl_matrix_complex_calloc( 2, 2 );

	gsl_matrix_complex_set ( Rot, 1, 1, gsl_complex_rect( cos(OMEGA*t/2), 0) ); 
	gsl_matrix_complex_set ( Rot, 1, 2, gsl_complex_rect( -sin(OMEGA*t/2), 0));
	gsl_matrix_complex_set ( Rot, 2, 1, gsl_complex_rect( sin(OMEGA*t/2), 0) );
	gsl_matrix_complex_set ( Rot, 2, 2, gsl_complex_rect( cos(OMEGA*t/2), 0) );

	/* inverse of Rot:
	 *
	 * first perform the LU decomposition of Rot: */	
	gsl_permutation* perm = gsl_permutation_calloc(2); 
	int signum;
	gsl_matrix_complex* Rot_lu = gsl_matrix_complex_calloc( 2, 2 );
	gsl_matrix_complex_memcpy ( Rot_lu, Rot );
	gsl_linalg_complex_LU_decomp ( Rot_lu, perm, &signum);
	/* then calculate Rot_inv: */
	gsl_matrix_complex* Rot_inv = gsl_matrix_complex_calloc( 2, 2 );
	gsl_linalg_complex_LU_invert ( Rot_lu, perm, Rot_inv );

	/* Finally calculate Rot * m * Rot_inv */
	gsl_blas_zgemm ( CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.,0),
			Rot, M, gsl_complex_rect(0,0), M );
	gsl_blas_zgemm ( CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.,0),
			M, Rot_inv, gsl_complex_rect(0,0), M );

	return (M);
}		/* -----  end of function rotate  ----- */


/* 
 *      FUNCTION  
 *         Name:  D_psi
 *  Description:  the matrix D[|psi><psi|]
 * 
 */
gsl_matrix_complex* D_psi ( gsl_matrix_complex* dpsi )
{
	gsl_matrix_complex* dp = gsl_matrix_complex_calloc(2, 2);

	return (dp);
}		/* -----  end of function D_psi  ----- */
