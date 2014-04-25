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
 * =====================================================================================
 *
 *       Filename:  Expi.c
 *
 *    Description:  Calculation of the exponential integrals (w/o the GSL special
 *    			functions). 
 *
 *        Version:  1.0
 *        Created:  15/04/2013 11:54:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"
#include <gsl/gsl_integration.h>
/* #include <gsl/gsl_sf_expint.h> */

double fu ( double t, void* params ) 
{
	double f = -exp(-t)/t ;
	return f ;
}

double ex ( double t, void* params )
{
	double e = exp(-t) ;
	return e ;
}

int expi ( double x, double* result, double* abserr )
{
	double r, err ;

	gsl_integration_workspace *expi_ws =
	       	gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_function F ;
	F.function = &fu ;

	int status ;
	status = gsl_integration_qagiu ( &F, x, 10e-7, .01 , WS_SZ, expi_ws, &r,
				&err) ;

	*result = r ; *abserr = err ;

/*	Using the GSL special functions, it is simply:
 *
 *      *result = - gsl_sf_expint_E1(x) ;
 */

	gsl_integration_workspace_free (expi_ws) ;

	return status;
}

int expi_plus ( double x, double* result, double* abserr )
{
	double r, err ;

	gsl_integration_workspace *expi_ws =
	       	gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_function F ;
	F.function = &ex ;

	int status ;
	status = gsl_integration_qawc ( &F, -x, x, 0, 1e-6, .01, WS_SZ,
			expi_ws, &r, &err ) ;

	double R , EXPI , ERREXPI ;
	int s = expi ( x , &EXPI, &ERREXPI ) ;
	R = - (r - EXPI) ;                          /* The minus (-) sign because of the def.
							of ex 				*/
	double ERR ;
	ERR = err + ERREXPI ;

	*result = R ; 
	*abserr = ERR ;

/* 	Using the GSL exponential integral, it is simply:
 *
 * 	R = gsl_sf_expint_Ei(x)  ;
 * 
 */
	gsl_integration_workspace_free (expi_ws) ;

	return status + s ;
}


