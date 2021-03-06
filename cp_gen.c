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
 *       Filename:  cp_gen.c
 *
 *    Description:  Bloch matrix for the completely positive generator
 *
 *        Version:  1.0
 *        Created:  28/04/2014 23:09:21
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "funcs.h"


/* 
 *      FUNCTION  
 *         Name:  cp_mat
 *  Description:  Creating the CP generator matrix in Bloch form, taking as
 *  		  arguments all the integrals, the physical parameters and 
 *  		  the pointer to the matrix
 * 
 */
int cp_mat ( gsl_matrix* cp_mx , void* params )
{
	double integrals[12] ;
	if ( (integration(params,integrals)) != 0 )
		return -1;

	/* Setting the integrals */
	double rcc, rss, rsc, isc, rcs, ics, rc0 ;
	rcc = integrals[0] ;
	rss = integrals[2] ;
	rsc = integrals[4] ;
	isc = integrals[5] ;
	rcs = integrals[6] ;
	ics = integrals[7] ;
	rc0 = integrals[8] ;

	/* Copying the parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	double D = sqrt(POW_2(o_1)-POW_2(O)) ;

	/* Ensuring that the matrix is zeroed */
	gsl_matrix_set_zero ( cp_mx ) ;

	/*  Building the matrix */
	unsigned int i ;
	for ( i = 0; i < 4; i++ )
		gsl_matrix_set ( cp_mx, 0 , i , 0 ) ;
	gsl_matrix_set ( cp_mx, 1, 0, 0 ) ; 
	gsl_matrix_set ( cp_mx, 1, 1, 2*(O/o_1)*rss+(1+POW_2(O/o_1))*rcc+2*POW_2(D/o_1)*rc0 ) ;
	gsl_matrix_set ( cp_mx, 1, 2, o_1/4-(1+POW_2(O/o_1))*rcs+2*(O/o_1)*rsc ) ;
	gsl_matrix_set ( cp_mx, 1, 3, 0	) ;
	gsl_matrix_set ( cp_mx, 2, 0, 0	) ;
	gsl_matrix_set ( cp_mx, 2, 1, -(o_1/4-(1+POW_2(O/o_1))*rcs+2*(O/o_1)*rsc) ) ;
	gsl_matrix_set ( cp_mx, 2, 2, 2*(O/o_1)*rss+(1+POW_2(O/o_1))*rcc+2*POW_2(D/o_1)*rc0 ) ;
	gsl_matrix_set ( cp_mx, 2, 3, 0	) ;
	gsl_matrix_set ( cp_mx, 3, 0, (1+POW_2(O/o_1))*ics-2*(O/o_1)*isc ) ;
	gsl_matrix_set ( cp_mx, 3, 1, 0	) ;
	gsl_matrix_set ( cp_mx, 3, 2, 0	) ;
	gsl_matrix_set ( cp_mx, 3, 3, (1+POW_2(O/o_1))*rcc+2*(O/o_1)*rss )  ;

	/* Multiplies times -4 to obtain the Bloch matrix */
	gsl_matrix_scale ( cp_mx, -4 ) ;

	return 0;
}		/* -----  end of function cp_mat  ----- */


