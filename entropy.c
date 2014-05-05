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

#include <math.h>
#include <gsl/gsl_sf_log.h>

#include "funcs.h"

/* 
 *      FUNCTION  
 *         Name:  entropy_production
 *  Description:  
 * 
 */
double entropy_production ( const gsl_vector* rho, const gsl_vector* rhoeq, const gsl_matrix* L  )
{
	/* l1, l2, l3 */
	double l[3] ; int i, j ;
	for ( i = 1 ; i < 3 ; i++ )
	{
		l[i] = 0 ;
		for ( j = 0 ; j < 3 ; j++ )
			l[i] += gsl_matrix_get(L,i,j)*gsl_vector_get(rho,j) ;
	}	

	/* L[rho] */
	double Lr = 0 ;
	for ( i = 1 ; i < 3 ; i++ )
		Lr += l[i]*gsl_vector_get(rho,i) ;

	/* L[rhoeq] */
	double Leq = 0 ;
	for ( i = 1 ; i < 3 ; i++ )
		Leq += l[i]*gsl_vector_get(rhoeq,i) ;

	/* r , req */
	double r, req ;
	r = req = 0 ;

	for ( i = 1 ; i < 3 ; i++ )
		r += POW_2(gsl_vector_get(rho,i)) ;
	r = sqrt(r) ;

	for ( i = 1 ; i < 3 ; i++ )
		req += POW_2(gsl_vector_get(rhoeq,i)) ;
	req = sqrt(req) ;

	/* internal entropy s */
	double s ;
	s = -(gsl_sf_log((1+r)/(1-r))*Lr/r - gsl_sf_log((1+req)/(1-req))*Leq/req) ;

	return s;
}		/* -----  end of function entropy_production  ----- */
