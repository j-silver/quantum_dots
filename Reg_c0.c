/*
 * =====================================================================================
 *
 *       Filename:  Reg_c0.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/04/2013 11:26:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), argenti@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int re_gc0 ( void* params, double* val)
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double vtemp = M_PI_2*alpha*(O*exp(-O/o_c)) ;
	double v = vtemp / tanh(b*O/2) ;

	*val = v ;
	
	return 0;
}
