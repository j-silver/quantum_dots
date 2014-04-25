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
 *       Filename:  Re_gss.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  14/04/2013 22:13:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int re_gss ( void* params, double* val )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;		
	alpha =  pars->alpha ;

	double v = M_PI_4*alpha*( (o_1-O)*exp(-(o_1-O)/o_c)/tanh(b*(o_1-O)/2)
				    - (o_1+O)*exp(-(o_1+O)/o_c)/tanh(b*(o_1+O)/2) ) ;
	*val = v ;
	
	return 0;
}
