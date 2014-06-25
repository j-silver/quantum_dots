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
 *       Filename:  sample.c
 *
 *    Description:  Sample random points on the sphere to check positivity 
 *    			violations
 *
 *        Version:  1.0
 *        Created:  16/05/2014 21:26:24
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include	<gsl/gsl_rng.h>
#include	<time.h>
#include	<math.h>
#include	<gsl/gsl_vector.h>
#include	<gsl/gsl_matrix.h>
#include	<stdlib.h>
#include	"funcs.h"

/* 
 *      FUNCTION  
 *         Name:  sample
 *  Description:  We choose N points in the Theta interval [0, Pi) and
 *  			N/sin(Theta) points in the Phi interval [0, 2Pi) .
 *  			
 *  			This is because the circles at colatitude Theta have
 *  			radius sin(Theta) instead of 1 and we must reduce the 
 *  			number of points accordingly, to have a uniform 
 *  			distribution all over the sphere.
 * 
 */
int sample ( const gsl_matrix* M, unsigned int N )
{
	if ( N == 0 )
		return -1;

	gsl_rng_env_setup() ;

	/* Double precision RANLXS generator */
	gsl_rng* r = gsl_rng_alloc(gsl_rng_ranlxd2) ;

	/* Taking the time in seconds as seed */
 	unsigned long int seed = (unsigned long int) time(0) ;

	/*  seeding  */
	gsl_rng_set(r, seed) ;

	/* Generating and testing points */
	double Theta, Phi ;
	gsl_vector* v = gsl_vector_calloc(4) ;
	double time_der ;

	FILE* f = fopen( "POS_VIOLATIONS", "w" );

	unsigned int i ;
	for ( i = 0; i < N; i++ )
	{
		Theta = M_PI*gsl_rng_uniform(r);
		Phi   = 2*M_PI*gsl_rng_uniform(r);
		bloch_vector(v, 1, Theta, Phi) ;
		time_der = r0_dot( M, v ) ;
		if ( time_der > 0 )
		{
			fprintf( f, "%-.9f %-.9f %-.9f %g\n",
					VECTOR(v,1),VECTOR(v,2),VECTOR(v,3),
					time_der ) ;
		}
	}
	
	fclose (f) ;

	gsl_rng_free(r) ;

	printf("n. of points: %d\n", N) ;

	return 0;
}		/* -----  end of function sample  ----- */



/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	unsigned int num = 1e6 ;
	
	gsl_matrix* red_m = gsl_matrix_calloc(4, 4) ;
	mat_read (red_m, "REDFIELD_MATRIX") ;

	sample (red_m, num) ;

	gsl_matrix_free(red_m) ;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
