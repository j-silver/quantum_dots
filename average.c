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
 *       Filename:  average.c
 *
 *    Description:  Integral of the entropy production for the Redfield dynamics
 *
 *        Version:  1.0
 *        Created:  31/05/2014 15:10:24
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include	<stdlib.h>
#include	<gsl/gsl_integration.h>
#include 	<gsl/gsl_errno.h>
#include	<string.h>

#include 	"funcs.h"
#include 	"initial.h"


/* 
 *      FUNCTION  
 *         Name:  entropy_prod_time
 *  Description:  Entropy production as a function of time t 
 * 
 */
double entropy_prod_time ( double t1, void* pars )
{
	/* Setting the parameters */
	struct ent_pars* p = (struct ent_pars*) pars ;
	const gsl_vector* rhoeq = p->rhoeq ;
	gsl_matrix* red_m = p->red_m ;

	/* Value to be returned */
	double ept ;

	/* 
	 *
	 * Initializing the system for Redfield dynamics
	 *
	 */
	gsl_odeiv2_system r_sys = { generator, jac, 4, (void*) red_m } ;

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	gsl_odeiv2_step* r_s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 ) ;

	/* Setting the step control: abserr=1e-9, relerr=1e-3 */
	gsl_odeiv2_control* r_c = gsl_odeiv2_control_standard_new ( 1e-9, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* r_e = gsl_odeiv2_evolve_alloc ( 4 ) ;

	/* Setting the initial state */
	gsl_vector* state = gsl_vector_calloc(4) ;
	unsigned int i ;                                 /* counter */
	for ( i = 0; i < 4; i++ )
		gsl_vector_set(state, i, R[i]) ;

	/* Evolving */
	double time = 0 ;
	while ( time < t1 )
	{
		evol( time, state, STEP, r_e, r_c, r_s, &r_sys ) ;
		time += STEP ;
	}
	ept = entropy_production( state, rhoeq, red_m ) ;

	/* free memory for evolution */
        gsl_odeiv2_evolve_free (r_e) ;
        gsl_odeiv2_control_free (r_c) ;
        gsl_odeiv2_step_free (r_s) ;

	return (ept) ;
}		/* -----  end of function entropy_prod_time  ----- */


/* 
 *      FUNCTION  
 *         Name:  entropy_integ
 *  Description:  Integrate entropy production up to time t
 * 
 */
double entropy_integ ( double t, const gsl_vector* rhoeq, gsl_matrix* red_m )
{
	double res, abserr;
	
	/*  Parameters to be passed to the integrand */	
	struct ent_pars params ;
	params.rhoeq = rhoeq ;
	params.red_m = red_m ;
	void* pars = (void*) &params ;

	/* Allocating workspace for integration */
	gsl_integration_workspace* ent_t = gsl_integration_workspace_alloc(1e6) ;

	/* Defining the integration function */
	gsl_function Ent ;
	Ent.function = &entropy_prod_time ;
	Ent.params   = pars ;

	/* Adaptive integration */
	gsl_integration_qag( &Ent, 0, t, 1e-3, 1e-1, 1e6, GSL_INTEG_GAUSS61,
			ent_t, &res, &abserr) ;

	return (res);
}		/* -----  end of function entropy_integ  ----- */


/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	gsl_matrix* red_m = gsl_matrix_calloc(4,4) ;
	mat_read ( red_m, "REDFIELD_MATRIX" ) ;

	gsl_vector* rhoeq = gsl_vector_calloc(4) ;
	FILE* f = fopen( "RED_STATIONARY.dat", "r") ;
	gsl_vector_fscanf ( f, rhoeq ) ;
	if ( f == NULL )
                printf("Error: %s.\nFailed to open RED_STATIONARY.dat\n",
				strerror(errno)) ;
	fclose (f) ;

	double INT = entropy_integ ( t_end, rhoeq, red_m ) ;
	printf("INTEGRATED ENT.: %g\n", INT) ;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
