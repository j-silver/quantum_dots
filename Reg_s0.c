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
 *       Filename:  Reg_s0.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/03/2014 21:10:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include <gsl/gsl_integration.h>

#include "funcs.h"

/* 
   We split the interval of integration (0,+inf) into three parts:
   1) (0, O/2)
   		and we use the diffeomorphism k -> 1/k to cope with the singularity
		in k = 0 ;
   2) (O/2, 3*O/2)
   		and we calculate the Cauchy principal value around k = O ;
   3) (3*O/2, +inf)
 */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fu_inv
 *  Description:  To integrate on the 1) interval we make the change of variable 
 *		  k -> 1/k 
 * =====================================================================================
 */
double fu_inv ( double k, void* params )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double temp = alpha*O*exp(-1/(k*o_c))/(1-(O*O*k*k)) ;
	double val = temp/(k*tanh(b/(k*2))) ;

	return val ;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fu_coth
 *  Description:  The function to be integrated around the singularity O, by calculating
 *		  the Cauchy principal value. We use qawc and therefore we multiply by
 *		  (k-O).
 * =====================================================================================
 */
double fu_cau ( double k, void* params )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double temp = alpha*O*k*exp(-k/o_c)/(k+O) ;
	double val = temp/tanh(b*k/2) ;

	return val ;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fu_coth
 *  Description:  Function to be integrated on the tail (3*O/2 , +inf)
 * =====================================================================================
 */
double fu_coth ( double k, void* params )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double temp = alpha*O*k*exp(-k/o_c)/(k*k-O*O) ;
	double val = temp/tanh(b*k/2) ;

	return val ;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  re_gs0
 *  Description:  Re g_s0
 * =====================================================================================
 */
int re_gs0 ( void* params, double* val, double* error )
{
	struct f_params* pars = (struct f_params*) params ;
	double O ;
	O = pars->Omega ;

	double r, err ;
	double r1, r2, r3, err1, err2, err3 ;

	gsl_function f, F, G ;
	f.function = &fu_coth ; F.function = &fu_cau ; G.function = &fu_inv ;
	f.params = pars ; F.params = pars ; G.params = pars ;

	gsl_integration_workspace *fu_fin_ws =
		gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_integration_workspace *fu_inv_ws =
		gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_integration_workspace *fu_inf_ws =
		gsl_integration_workspace_alloc (WS_SZ) ;

	int status1 ;
	status1 = gsl_integration_qagiu ( &G, 2/O, 1e-9, .001, WS_SZ,
		       	fu_inv_ws, &r1,	&err1 ) ;
	int status2 ;
	status2 = gsl_integration_qawc ( &F, O/2, 3*O/2, O, 1e-9, .001, WS_SZ,
			fu_fin_ws, &r2, &err2 ) ;

	int status3 ;
	status3 = gsl_integration_qagiu ( &f, 3*O/2, 1e-9, .001, WS_SZ,
		       	fu_inf_ws, &r3,	&err3 ) ;

	r = r1 + r2 + r3 ; err = err1 + err2 + err3 ;
	*val = r ; *error = err ;

	gsl_integration_workspace_free(fu_fin_ws) ;
	gsl_integration_workspace_free(fu_inv_ws) ;
	gsl_integration_workspace_free(fu_inf_ws) ;

	int status = status1 + status2 + status3 ;
	return status ;
}


