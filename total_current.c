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
 *       Filename:  total_current.c
 *
 *    Description:  Total current ( AC + DC )
 *
 *        Version:  1.0
 *        Created:  26/06/2014 14:59:12
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "funcs.h"
#include <gsl/gsl_vector.h>

/* 
 *      FUNCTION  
 *         Name:  tot_current
 *  Description:  
 * 
 */
double tot_current ( const gsl_vector* rho, void* params )
{
	struct f_params* pars = (struct f_params*) params ;
	double Omega = pars->Omega ;
	double o_1   = pars->omega_1 ;
	double D     = sqrt(o_1*o_1 - Omega*Omega) ;
	double curr = -((VECTOR(rho, 2)*D + VECTOR(rho, 3)*Omega)/o_1) ;
	return (curr) ;
}		/* -----  end of function tot_current  ----- */
