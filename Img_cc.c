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
 *       Filename:  Img_cc.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/04/2013 17:08:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), argenti@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int im_gcc ( void* params, double* val, double* error )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double e1, e2, E1, E2 ;
	double err1, err2, ERR1, ERR2 ;

	double beta1 = O + o_1 ;
	double beta2 = o_1 - O ;
	double mu = 1/o_c ;

	expi ( beta1*mu, &e1, &err1 ) ;
	expi_plus ( beta1*mu, &E1, &ERR1 ) ;
	expi ( beta2*mu, &e2, &err2 ) ;
	expi_plus ( beta2*mu, &E2, &ERR2 ) ;

	double imgcc = -alpha/mu - alpha/4 *
			( + beta1*exp(beta1*mu)*e1 + beta2*exp(beta2*mu)*e2
			- beta1*exp(-beta1*mu)*E1 - beta2*exp(-beta2*mu)*E2 ) ;
  	*val = imgcc ; 

	double err =  alpha/4 *
			( + beta1*exp(beta1*mu)*err1 + beta2*exp(beta2*mu)*err2
			+ beta1*exp(-beta1*mu)*ERR1 + beta2*exp(-beta2*mu)*ERR2 ) ;

	*error = err ;

	return 0;
}
