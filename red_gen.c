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
 *       Filename:  red_gen.c
 *
 *    Description:  Redfiled-type Generator 
 *
 *        Version:  1.0
 *        Created:  25/04/2014 22:28:19
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
 *         Name:  red_mat
 *  Description:  Creating the Redfield generator matrix in Bloch form, taking as
 *  		  arguments all the integrals, the physical parameters 
 *  		  and the pointer to the matrix
 * 
 */
int red_mat ( gsl_matrix* red_mx , void* params )
{
	double integrals[12] ;
	if ( integration(params,integrals) != 0 )
		return -1;

	/* Setting the integrals */
	double rcc, icc, rss, iss, rsc, isc, rcs, ics, rc0, ic0, rs0, is0 ;
	rcc = integrals[0] ;
	icc = integrals[1] ;
	rss = integrals[2] ;
	iss = integrals[3] ;
	rsc = integrals[4] ;
	isc = integrals[5] ;
	rcs = integrals[6] ;
	ics = integrals[7] ;
	rc0 = integrals[8] ;
	ic0 = integrals[9] ;
	rs0 = integrals[10] ;
	is0 = integrals[11] ;

	/* Copying the parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	double D = sqrt(POW_2(o_1)-POW_2(O)) ;

	/* Ensuring that the matrix is zeroed */
	gsl_matrix_set_zero ( red_mx ) ;

	/* Building the matrix */
	gsl_matrix_set ( red_mx, 1, 0, (D/o_1)*(-(O/o_1)*ic0+iss+(O/o_1)*icc) ) ;
	gsl_matrix_set ( red_mx, 1, 1, POW_2(D/o_1)*rc0+(O/o_1)*rss+POW_2(O/o_1)*rcc ) ;
	gsl_matrix_set ( red_mx, 1, 2, o_1/4+(O/o_1)*rsc-POW_2(O/o_1)*rcs ) ;
	gsl_matrix_set ( red_mx, 1, 3, (D/o_1)*(rsc-(O/o_1)*rcs) ) ;
	gsl_matrix_set ( red_mx, 2, 0, (D/o_1)*(is0+isc-(O/o_1)*ics) ) ;
	gsl_matrix_set ( red_mx, 2, 1, -(o_1/4)-(O/o_1)*rsc+rcs ) ;
	gsl_matrix_set ( red_mx, 2, 2, POW_2(D/o_1)*rc0+rcc+(O/o_1)*rss ) ;
	gsl_matrix_set ( red_mx, 2, 3, -(D/o_1)*((O/o_1)*rcc+rss) ) ;
	gsl_matrix_set ( red_mx, 3, 0, (1+POW_2(O/o_1))*ics-2*(O/o_1)*isc ) ;
	gsl_matrix_set ( red_mx, 3, 1, -(D/o_1)*rs0 ) ; 
	gsl_matrix_set ( red_mx, 3, 2, -(O*D/(o_1*o_1))*rc0 ) ;
	gsl_matrix_set ( red_mx, 3, 3, (1+POW_2(O/o_1))*rcc+2*(O/o_1)*rss ) ;

	/* Multiplies times -2 to obtain the Bloch matrix */
	gsl_matrix_scale ( red_mx, -2 ) ;


	return 0;
}		/* -----  end of function red_mat  ----- */


