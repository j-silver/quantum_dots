/*
 * =====================================================================================
 *
 *       Filename:  cos_coth.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  20/04/2013 15:27:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri
 *   Organization:  
 *
 * =====================================================================================
 */

#include <gsl/gsl_integration.h>

#include "funcs.h"
#include <stdio.h>

double k_integrand ( void* params, double k, double w )
{
	struct ext_pars* pars = (struct ext_pars*) params ;

	double o_c, b, alpha ;
	o_c   = (pars->params.omega_c) ;
	b     = (pars->params.beta) ;
	alpha = (pars->params.alpha) ;

	/* break into 2 parts, trying to avoid round-off induced errors */
	double temp = alpha*k*exp(-k/o_c)*cos(k*w) ;
	double val  = temp/tanh(b*k/2) ;
	return val ;
}

double k_func ( double k, void* x_ps )
{
	struct ext_pars* pr = (struct ext_pars*) x_ps ;
	
	double K = k_integrand(pr, k, pr->w) ;
	return K ;
}

double k_integration ( double w, void* params )
{
	struct ext_pars ex_p ;
	struct f_params* p = (struct f_params *) params ;
	ex_p.params = *p ;
	ex_p.w = w ;

	gsl_function K ;
	K.function = &k_func ;
	K.params = &ex_p ;

	/* allocating memory in workspace (limit = 1,000,000) */
	gsl_integration_workspace* k_alloc = 
		gsl_integration_workspace_alloc(1E+6) ;

	double result, abserr ;
	gsl_integration_qagiu ( &K, 0, 1e-6, 1e-2, 1E+6, k_alloc, &result,
				&abserr ) ;

	gsl_integration_workspace_free (k_alloc) ;

	return result;
}

double w_integrand ( double w, void* params )
{
	struct f_params* pars = (struct f_params*) params ;
	double O, o_1 ;
	O     = (pars->Omega) ;
	o_1   = (pars->omega_1) ;

	double integ = k_integration(w, pars) ;
	double W = sin(-O*w)*cos(-o_1*w)*integ ;
	return W ;
}
	
int w_integration ( double* result, double* abserr, void* params )
{
	struct f_params* pars = (struct f_params *) params ;

	gsl_function W ;
	W.function = &w_integrand ;
	W.params = pars ;

	/* allocating memory in workspace */
	gsl_integration_workspace* w_alloc =
		gsl_integration_workspace_alloc(1e6) ;

	double r, err ;
	gsl_integration_qag ( &W, 0, 10, 1e-6, 1e-2, 1e6, GSL_INTEG_GAUSS61,
			w_alloc, &r, &err ) ;	

	gsl_integration_workspace_free (w_alloc) ;

	*result = r ;
	*abserr = err ;

	return 0 ;
}

