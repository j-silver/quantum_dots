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
 *       Filename:  polar.c
 *
 *    Description:  Cartesian <-> Polar coordinates transformation
 *
 *        Version:  1.0
 *        Created:  06/05/2014 14:21:40
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

/* 
 *      FUNCTION  
 *         Name:  bloch_vector
 *  Description:  Return the Bloch vector corresponding to (r,theta,phi) polar coords
 * 
 */
int bloch_vector ( gsl_vector* v, unsigned double r, double theta, double phi )
{
	/* Check the validity of coordinates */
	if ( r > 1 )                   /* Non physical state */
	{
		v = NULL ;
		return -1;
	}
	theta = fmod (theta, M_PI ) ;          /* Taking the modulo */
	phi   = fmod ( phi, 2*M_PI )  ;

	gsl_vector_set ( v, 1, r*sin(theta)*cos(phi) ) ;
	gsl_vector_set ( v, 2, r*sin(theta)*sin(phi) ) ;
	gsl_vector_set ( v, 3, r*cos(theta) ) ;

	return 0;
}		/* -----  end of function bloch_vector  ----- */


/* 
 *      FUNCTION  
 *         Name:  polars
 *  Description:  From the Bloch vector representation to the polar coordinates
 * 
 */
int polars ( unsigned double* r, double* theta, double* phi, const gsl_vector* v )
{
	*r = gsl_hypot3( gsl_vector_get(v,1), gsl_vector_get(v,2), gsl_vector_get(v,3) ) ;
	if ( *r > 1 )
		return -1 ;

	*theta = acos(gsl_vector_get(v,3)) ;

	if ( fabs(*theta)==1 )
		*phi = 0 ;
	else
		*phi = atan2(gsl_vector_get(v,2),gsl_vector_get(v,1)) ;

	return 0;
}		/* -----  end of function polars  ----- */
