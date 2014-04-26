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
 *       Filename:  integs.c
 *
 *    Description:  Populating the integrals array
 *
 *        Version:  1.0
 *        Created:  26/04/2014 13:56:01
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "funcs.h"

/* 
 *      FUNCTION  
 *         Name:  integrals
 *  Description:  Populate the integrals array 
 * 
 */
int integration ( void* params, double* integrals )
{
/* 	struct f_params* pars = (struct f_params*) params ;
 *         double o_c, b, O, o_1, alpha ;
 *         assign_p ( pars, &o_c, &b, &O, &o_1) ;
 *         alpha = pars->alpha ;	
 * 
 */
	double regcc ;
	re_gcc( params, &regcc ) ;
	integrals[0] = regcc ;

	double imgcc, imgcc_error ;
	im_gcc ( params, &imgcc, &imgcc_error ) ;
	integrals[1] = imgcc ;

	double regss ;
	re_gss ( params, &regss ) ;	
	integrals[2] = regss ;

	double imgss, imgss_error ;
	im_gss ( params, &imgss, &imgss_error ) ;
	integrals[3] = imgss ;

	double regsc, regscerr ;
	re_gsc ( &regsc, &regscerr, params ) ;
	integrals[4] = regsc ;

	double imgsc ;
	im_gsc ( params, &imgsc ) ;
	integrals[5] = imgsc ;

	double regcs, regcserr ;
	re_gcs ( &regcs, &regcserr, params ) ;
	integrals[6] = regcs ;

	double imgcs ;
	im_gcs ( params, &imgcs ) ;
	integrals[7] = imgcs ;

	double regc0 ;
	re_gc0 ( params, &regc0 ) ;
	integrals[8] = regc0 ;

	double imgc0, imgc0_error ;
	im_gc0 ( params, &imgc0, &imgc0_error ) ;
	integrals[9] = imgc0 ;

	double regs0, regs0_error ;
	re_gs0 ( params, &regs0, &regs0_error ) ;
	integrals[10] = regs0 ;

	double imgs0 ;
	im_gs0 ( params, &imgs0 ) ;
	integrals[11] = imgs0 ;

	return 0;
}		/* -----  end of function integrals  ----- */
