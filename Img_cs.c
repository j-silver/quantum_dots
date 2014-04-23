/*
 * =====================================================================================
 *
 *       Filename:  Im_gcs.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  14/04/2013 22:49:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author: Giuseppe Argentieri
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int im_gcs ( void* params, double* val )
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;		
	alpha =  pars->alpha ;

	double v = M_PI_4*alpha*((o_1-O)*exp(-(o_1+O))/o_c 
		+ (o_1-O)*exp(-(o_1-O)/o_c)) ;

	*val = v ;

	return 0;
}
