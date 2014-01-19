/*
 * =====================================================================================
 *
 *       Filename:  Img_c0.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/04/2013 13:43:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), argenti@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"
#include <stdio.h>

int im_gc0 ( void* params, double* val, double* error ) 
{
	struct f_params* pars = (struct f_params*) params ;
	double o_c, b, O, o_1, alpha ;
	assign_p ( pars, &o_c, &b, &O, &o_1 ) ;
	alpha = pars->alpha ;

	double ei, Ei, err1, err2 ;
	expi( O/o_c, &ei, &err1 ) ;
	expi_plus ( O/o_c, &Ei, &err2 ) ;

	double v = -alpha*M_PI*(2*o_c-O/2*(exp(O/o_c)*(-ei) - exp(-O/o_c)*(Ei+ei))) ;
	*val = v ;
	double err = (alpha*M_PI*O/2) * (exp(O/o_c)*err1 - exp(-O/o_c)*err2) ; 
	*error = err ;

	return 0;
}
