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
 *       Filename:  red_evol.c
 *
 *    Description:  Evolution of the Redfield dynamics
 *
 *        Version:  1.0
 *        Created:  15/05/2014 19:14:06
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include	<stdlib.h>
#include	<string.h>
#include	<gsl/gsl_ieee_utils.h>
#include 	"funcs.h"
#include	"initial.h"


/* 
 *      FUNCTION  
 *         Name:  red_evol
 *  Description:  
 * 
 */
int red_evol ( void* params, const double r[], double time_end, double step,
		const gsl_vector* req_red, gsl_matrix* red_m )
{
	struct f_params* pars = (struct f_params*) params ;
	double Omega, omega_1 ;
	Omega   = pars->Omega ;
	omega_1 = pars->omega_1 ;

	unsigned int i ;                        /* counter for the for loops */

	/* 
	 *
	 *  evolving the systems
	 * 
	 */
	double t = 0 ; 

	/* setting the initial vector */
	gsl_vector* init_red = gsl_vector_calloc(4) ;
	for ( i = 0 ; i < 4 ; i++ )
		gsl_vector_set ( init_red, i, r[i] ) ;

	/* 
	 *
	 * Initializing the system for Redfield dynamics
	 *
	 */
	gsl_odeiv2_system red_sys = { generator, jac, 4, (void*) red_m } ;

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	gsl_odeiv2_step* r_s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 ) ;

	/* Setting the step control: abserr=1e-9, relerr=1e-3 */
	gsl_odeiv2_control* r_c = gsl_odeiv2_control_standard_new ( 1e-9, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* r_e = gsl_odeiv2_evolve_alloc ( 4 ) ;

	/* opening the files */
	FILE* f_red = fopen ( "RED-EVOLUTION.dat", "w" ) ;
	FILE* g_red = fopen ( "RED-ENTROPY.dat", "w" ) ;
	FILE* h_red = fopen ( "RED-CURRENT.dat", "w" ) ;


	/* writing data */
	while ( t < t_end )
	{
		evol ( t, init_red, step, r_e, r_c, r_s, &red_sys ) ;
		fprintf ( f_red, "%.2f %.9f %.9f %.9f %.9f\n", t,
				VECTOR(init_red,1), VECTOR(init_red,2),
				VECTOR(init_red,3),
				gsl_hypot3(VECTOR(init_red,1),VECTOR(init_red,2),
					VECTOR(init_red,3)) ) ;
		fprintf ( g_red, "%.2f %.9f\n",
				t, entropy_production( init_red, req_red, red_m )) ;
		fprintf ( h_red, "%.2f %.9f\n",
				t, -VECTOR(init_red,3)*Omega/omega_1 ) ;
		t += step ;
	}
	
	/* final entropy */
	printf("Final entropy: %g\n", entropy_of_state(init_red)) ;

	/*  close the files */
	fclose (f_red) ;
	fclose (g_red) ;
	fclose (h_red) ;

	/* free memory for evolution */
	gsl_odeiv2_evolve_free (r_e) ;
	gsl_odeiv2_control_free (r_c) ;
	gsl_odeiv2_step_free (r_s) ;

	return 0;
}		/* -----  end of function red_evol  ----- */


/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	gsl_ieee_env_setup () ;			/* read GSL_IEEE_MODE */

	double beta = 1.0/T ;                   /* Boltzmann factor: beta */
	double omega_1 = gsl_hypot(OMEGA,D) ;   /* omega' */

	struct f_params params;
	params.omega_c = omega_c ;
	params.beta = beta ;
	params.Omega = OMEGA ;
	params.omega_1 = omega_1 ;
	params.alpha = alpha ;

	/* read the Redfield matrix from a file */
	gsl_matrix* red_m = gsl_matrix_calloc ( 4, 4 ) ;
	mat_read ( red_m, "REDFIELD_MATRIX" ) ;

	/* 
	 * Find the stationary state by solving the linear systems
	 * and the associated stationary currents
	 */
	gsl_vector* req_red = gsl_vector_calloc (4) ;
	FILE* f_red = fopen ( "RED_STATIONARY.dat", "r" ) ;
	if ( f_red == NULL )
		printf("error: %s\n", strerror(errno)) ;
	gsl_vector_fscanf ( f_red, req_red ) ;
	fclose ( f_red ) ;

	printf("REDFIELD DYNAMICS\n") ;
	printf("Stationary state: ( %.1f , %.9f , %.9f , %.9f )\n", 
			VECTOR(req_red,0), VECTOR(req_red,1),
			VECTOR(req_red,2), VECTOR(req_red,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.9f\n\n",
			-VECTOR(req_red,3)*OMEGA/omega_1) ;
	
	int status = red_evol ( &params, R, t_end, STEP, req_red, red_m ) ;
	if ( status != 0 )
		exit (EXIT_FAILURE) ;

	/* free memory for matrix */
	gsl_matrix_free(red_m) ;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */


