/*
 * =====================================================================================
 *
 *       Filename:  Img_ss.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/04/2013 16:44:51
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), argenti@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"

int im_gss ( void* params, double* val, double* error )
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
	expi_plus ( beta2*mu, &E2, &ERR2 ) ;
	expi ( beta2*mu, &e2, &err2 ) ;

	double imgss = (alpha/4)*(beta1*(e1*exp(beta1*mu))-E1*exp(-beta1*mu) +
				  beta2*(E2*exp(-beta2*mu)-e2*exp(beta2*mu))) ;
	*val = imgss ;

	double err = (alpha/4)*(beta1*(err1*exp(beta1*mu))-ERR1*exp(-beta1*mu) +
				  beta2*(ERR2*exp(-beta2*mu)-err2*exp(beta2*mu))) ;
	*error = err ;	

	return 0;
}



