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
 *       Filename:  Reg_sc.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/04/2014 12:27:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * =====================================================================================
 */

/*
   Re g_sc is made of 4 integrands: 2 with principal values and 2 without.

   The common factor function is k_func = alpha*k*exp(-k/omega_c)/tanh(b*k/2) 

   The poles are in (omega1-O) and (omega1+O) 

   The sign of the four fractions

   1/(k+(omega1+O)) , 1/(k+(omega1-O)) , 1/(k-(omega1-O)) , 1/(k-(omega1+O))

   are

   + , - , + , -
*/

#include <stdio.h>

#include "funcs.h"
#include <gsl/gsl_integration.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  k_func
 *  Description:  Common integrand factor
 * =====================================================================================
 */
double k_func ( double k, void* params )
{

	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double temp = alpha*k*exp(-k/o_c) ;
	double val = temp/tanh(b*k/2) ;
	
	return val;
}		/* -----  end of function k_func  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  k_func_1
 *  Description:  First integrand: k_func/(k+(omega1+O))
 * =====================================================================================
 */
double k_func_1 ( double k, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	/* Invoking k_func and returning value */
	double v , v1 ;
	v = k_func( k, pars ) ;
	v1 = v/(k+(o_1+O)) ;

	return v1 ;
}		/* -----  end of function k_func_1  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  first_int
 *  Description:  Integrating k_func/(k+(omega1+O))
 * =====================================================================================
 */
int first_int ( double* val, double* error, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;
	
	/* Allocating workspace for integration */
	gsl_integration_workspace* first_int_ws =
		gsl_integration_workspace_alloc(WS_SZ) ;

	gsl_function F ;
	F.function = &k_func_1 ;
	F.params = pars ;

	/* Integration over (0,+Infinity) */
	double res, err ;
	int status = gsl_integration_qagiu( &F, 0, 1e-9, 1e-3, WS_SZ, first_int_ws, 
			&res, &err ) ;

	*val = res ; *error = err ;

	gsl_integration_workspace_free(first_int_ws) ;
	return status;
}		/* -----  end of function first_int  ----- */


/* 
 *      FUNCTION  
 *         Name:  k_func_2
 *  Description:  Second integrand: k_func/(k+omega1-O)
 * 
 */
double k_func_2 ( double k , void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;
	
	/* Invoking k_func and returning value */
	double v, v2 ;
	v = k_func( k, pars ) ;
	v2 = v/(k+(o_1-O)) ;

	return v2;
}		/* -----  end of function k_func_2  ----- */


/* 
 *      FUNCTION  
 *         Name:  second_int
 *  Description:  Integration k_func/(k+(omega1-O))
 * 
 */
int second_int ( double* val, double* error, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	/* Allocating workspace for integration */
	gsl_integration_workspace* second_int_ws =
		gsl_integration_workspace_alloc(WS_SZ) ;

	gsl_function F;
	F.function = &k_func_2 ;
	F.params = pars ;

	double res, err ;
	int status = gsl_integration_qagiu( &F, 0, 1e-9, 1e-3, WS_SZ, second_int_ws,
			&res, &err ) ;
	*val = res ; *error = err ;

	gsl_integration_workspace_free(second_int_ws) ;

	return status;
}		/* -----  end of function second_int  ----- */


/* 
 *      FUNCTION  
 *         Name:  k_func_3
 *  Description:  k_func/(k-(omega1-O)) to be integrated on the tail +inf
 * 
 */
double k_func_3 ( double k, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	double temp , val ;
	temp = k_func( k, pars ) ;
	val = temp/(k-(o_1-O)) ; 

	return val;
}		/* -----  end of function k_func_3  ----- */


/* 
 *      FUNCTION  
 *         Name:  third_int
 *  Description:  Integration of the p.v. k_func/(k-(omega1-O))
 * 
 */
int third_int ( double* val, double* err, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	/* Integration of the p.v. */
	gsl_function F ;
	F.function = &k_func ;
	F.params = pars ;

	gsl_integration_workspace* third_int_ws =
		gsl_integration_workspace_alloc(WS_SZ) ;

	double v1 , err1 ;
	int status1 = gsl_integration_qawc( &F, (o_1-O)/2, 3*(o_1-O)/2, o_1-O, 
			1e-9, 1e-3, WS_SZ, third_int_ws, &v1, &err1 ) ;

	gsl_integration_workspace_free(third_int_ws) ;

	/* Integration to +infinity */
	gsl_function G ;
	G.function = &k_func_3 ;
	G.params = pars ;

	gsl_integration_workspace* third_int_ws_2 =
		gsl_integration_workspace_alloc(WS_SZ) ;

	double v2, err2 ;
	int status2 = gsl_integration_qagiu( &G, 3*(o_1-O)/2, 1e-9, 1e-3, WS_SZ,
			third_int_ws_2, &v2, &err2 ) ;

	gsl_integration_workspace_free(third_int_ws_2) ;

	/* Integration on ( 0, (o_1-O)/2 ) */
	gsl_integration_workspace* third_int_ws_3 = 
		gsl_integration_workspace_alloc(WS_SZ) ;

	gsl_function H ;
	H.function = &k_func_3 ;
	H.params = pars ;

	double v3, err3 ;
	int status3 = gsl_integration_qag( &H, 0, (o_1-O)/2, 1e-9, 1e-3, WS_SZ, 6,
			third_int_ws_3, &v3, &err3 ) ;

	gsl_integration_workspace_free(third_int_ws_3) ;

	*val = v1 + v2 + v3 ; *err = err1 + err2 + err3 ;

	return status1 + status2 + status3 ;
}		/* -----  end of function third_int  ----- */


/* 
 *      FUNCTION  
 *         Name:  k_func_4
 *  Description:  k_func/(k-(omega1+O))
 * 
 */
double k_func_4 ( double k, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;
	
	double temp , val ;
	temp = k_func( k, pars ) ;
	val = temp/(k-(o_1+O)) ; 

	return val;
}		/* -----  end of function k_func_4  ----- */


/* 
 *      FUNCTION  
 *         Name:  fourth_integ
 *  Description:  Integration of the p.v. kfunc/(k-(omega1+O))
 * 
 */
int fourth_int ( double* val, double* err, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	/* Integration of the p.v. */
	gsl_function F ;
	F.function = &k_func ;
	F.params = pars ;

	gsl_integration_workspace* fourth_int_ws =
		gsl_integration_workspace_alloc(WS_SZ) ;

	double v1 , err1 ;
	int status1 = gsl_integration_qawc( &F, (o_1+O)/2, 3*(o_1+O)/2, o_1+O, 
			1e-9, 1e-3, WS_SZ, fourth_int_ws, &v1, &err1 ) ;

	gsl_integration_workspace_free(fourth_int_ws) ;

	/* Integration to +infinity */
	gsl_function G ;
	G.function = &k_func_4 ;
	G.params = pars ;

	gsl_integration_workspace* fourth_int_ws_2 =
		gsl_integration_workspace_alloc(WS_SZ) ;

	double v2, err2 ;
	int status2 = gsl_integration_qagiu( &G, 3*(o_1+O)/2, 1e-9, 1e-3, WS_SZ,
			fourth_int_ws_2, &v2, &err2 ) ;

	gsl_integration_workspace_free(fourth_int_ws_2) ;

	/* Integration on ( 0 , (O+o_1)/2 ) */
	gsl_integration_workspace* fourth_int_ws_3 =
		gsl_integration_workspace_alloc(WS_SZ) ;

	gsl_function K ;
	K.function = &k_func_4 ;
	K.params = pars ;
	double v3, err3 ;

	int status3 = gsl_integration_qag( &K, 0, (O+o_1)/2, 1e-9, 1e-3, WS_SZ, 
			GSL_INTEG_GAUSS61, fourth_int_ws_3, &v3, &err3 ) ;

	gsl_integration_workspace_free(fourth_int_ws_3) ;

	*val = v1 + v2 + v3 ; *err = err1 + err2 + err3 ;

	return status1 + status2 + status3;
}		/* -----  end of function fourth_integ  ----- */


/* 
 *      FUNCTION  
 *         Name:  re_gsc
 *  Description:  
 * 
 */
int re_gsc ( double* result, double* error, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	double v1, v2, v3, v4, err1, err2, err3, err4 ;
	int status = first_int( &v1, &err1, pars ) + second_int( &v2, &err2, pars )
		+ third_int( &v3, &err3, pars ) + fourth_int( &v4, &err4, pars ) ;
	*result = -(v1 - v2 + v3 - v4)/4.0 ;
	*error = (err1 + err2 + err3 + err4)/4.0 ;

	return status;
}		/* -----  end of function re_gsc  ----- */

/* 
 *      FUNCTION  
 *         Name:  re_gcs
 *  Description:  Re g_cs is the same of Re g_sc but with O and o_1 exchanged. 
 *  		  Therefore the signs are
 *
 *  		  + , + , - , -
 * 
 */
int re_gcs ( double* result, double* error, void* params )
{
	/* Copying parameters */
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1 ;
	assign_p( pars, &o_c, &b, &O, &o_1 ) ;

	double v1, v2, v3, v4, err1, err2, err3, err4 ;
	int status =  first_int( &v1, &err1, pars ) + second_int( &v2, &err2, pars )
		+ third_int( &v3, &err3, pars ) + fourth_int( &v4, &err4, pars ) ;

	*result = -(v1 + v2 - v3 - v4)/4 ; *error = (err1 + err2 + err3 + err4)/4 ;

	return status;
}		/* -----  end of function re_gcs  ----- */
