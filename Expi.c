/*
 * =====================================================================================
 *
 *       Filename:  Expi.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/04/2013 11:54:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Giuseppe Argentieri (), argenti@ts.infn.it
 *   Organization:  
 *
 * =====================================================================================
 */

#include "funcs.h"
#include <gsl/gsl_integration.h>

int expi ( double x, double* result, double* abserr )
{
	double r, err ;

	gsl_integration_workspace *expi_ws =
	       	gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_function F ;
	F.function = &fu ;

	int status ;
	status = gsl_integration_qagiu ( &F, x, 10e-7, .01 , WS_SZ, expi_ws, &r,
				&err) ;

	*result = r ; *abserr = err ;
	
	gsl_integration_workspace_free (expi_ws) ;

	return status;
}

int expi_plus ( double x, double* result, double* abserr )
{
	double r, err ;

	gsl_integration_workspace *expi_ws =
	       	gsl_integration_workspace_alloc (WS_SZ) ;

	gsl_function F ;
	F.function = &ex ;

	int status ;
	status = gsl_integration_qawc ( &F, -x, x, 0, 1e-6, .01, WS_SZ,
			expi_ws, &r, &err ) ;

	*result = r ; 
	*abserr = err ;

	return status;
}

double fu ( double t, void* params ) 
{
	double f = exp(-t)/t ;
	return f ;
}

double ex ( double t, void* params )
{
	double e = exp(-t) ;
	return e ;
}
