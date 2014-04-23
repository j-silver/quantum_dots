/*
 * =====================================================================================
 *
 *       Filename:  Img_s0.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  19/04/2014 15:37:21
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), giuseppe.argentieri@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int im_gs0 ( void* params, double* val )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double temp = M_PI_2*alpha*O*exp(-O/o_c) ;

	*val = temp ;

	return 0;
}
