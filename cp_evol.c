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
 *       Filename:  cp_evol.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15/05/2014 19:36:38
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include	<gsl/gsl_odeiv2.h>
#include	<stdlib.h>
#include	<string.h>
#include	<gsl/gsl_ieee_utils.h>
#include	"funcs.h"
#include	"initial.h"

/* 
 *      FUNCTION  
 *         Name:  cp_evol
 *  Description:  
 * 
 */
int cp_evol ( void* params, const double r[], double time_end, double step,
		const gsl_vector* req_cp, gsl_matrix* cp_m )
{
	struct f_params* pars = (struct f_params*) params ;

	unsigned int i ;                        /* counter for the for loops */

	/* 
	 *
	 *  evolving the systems
	 * 
	 */
	double t = 0 ; 

	/* setting the initial vector */
	gsl_vector* init_cp = gsl_vector_calloc(4) ;
	for ( i = 0 ; i < 4 ; i++ )
		gsl_vector_set ( init_cp, i, r[i] ) ;

	/* 
	 *
	 * Initializing the system for Redfield dynamics
	 *
	 */
	gsl_odeiv2_system cp_sys = { generator, jac, 4, (void*) cp_m } ;

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	gsl_odeiv2_step* c_s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 ) ;

	/* Setting the step control: abserr=1e-9, relerr=1e-3 */
	gsl_odeiv2_control* c_c = gsl_odeiv2_control_standard_new ( 1e-9, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* c_e = gsl_odeiv2_evolve_alloc ( 4 ) ;

	/* opening the files */
	FILE* f_cp = fopen ( "CP-EVOLUTION.dat", "w" ) ;
	FILE* g_cp = fopen ( "CP-ENTROPY-PROD.dat", "w" ) ;
	FILE* h_cp = fopen ( "CP-CURRENT.dat", "w" ) ;
	FILE* i_cp = fopen ( "CP-ENTROPY.dat", "w" ) ;

	/* writing data */
	while ( t < t_end )
	{
		evol ( t, init_cp, step, c_e, c_c, c_s, &cp_sys ) ;
		fprintf ( f_cp, "%.2f %.9f %.9f %.9f %.9f\n", t,
				VECTOR(init_cp,1), VECTOR(init_cp,2),
				VECTOR(init_cp,3),
				gsl_hypot3(VECTOR(init_cp,1),VECTOR(init_cp,2),
					VECTOR(init_cp,3)) ) ;
		fprintf ( g_cp, "%.2f %.9f\n",
				t, entropy_production( init_cp, req_cp, cp_m )) ;
		fprintf ( h_cp, "%.2f %.9f\n",
				t, tot_current(init_cp, pars) ) ;
		fprintf ( i_cp, "%.2f %.9f\n",
				t, entropy_of_state(init_cp) ) ;
		t += step ;
	}

	/* final entropy */
	printf("Final entropy: %g\n", entropy_of_state(init_cp)) ;

	/*  close the files */
	fclose (f_cp) ;
	fclose (g_cp) ;
	fclose (h_cp) ;
	fclose (i_cp) ;

	/* free memory for evolution */
	gsl_odeiv2_evolve_free (c_e) ;
	gsl_odeiv2_control_free (c_c) ;
	gsl_odeiv2_step_free (c_s) ;

	return 0;
}		/* -----  end of function cp_evol  ----- */



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

	/* read the CP matrix from a file */
	gsl_matrix* cp_m = gsl_matrix_calloc ( 4, 4 ) ;
	mat_read ( cp_m, "CP_MATRIX" ) ;

	/* 
	 * Find the stationary state by solving the linear systems
	 * and the associated stationary currents
	 */

	gsl_vector* req_cp = gsl_vector_calloc (4) ;
	
	FILE* f_cp = fopen ( "CP_STATIONARY.dat", "r" ) ;
	gsl_vector_fscanf ( f_cp, req_cp ) ;
	if ( f_cp == NULL )
		printf("Error: %s.\nFailed to open CP_STATIONARY.dat",
				strerror(errno)) ;
	fclose ( f_cp ) ;

	printf("CP DYNAMICS\n") ;
	printf("Stationary state: ( %.1f , %.9f , %.9f , %.9f )\n", 
			VECTOR(req_cp,0), VECTOR(req_cp,1),
			VECTOR(req_cp,2), VECTOR(req_cp,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.9f\n\n",
			-VECTOR(req_cp,3)*OMEGA/omega_1) ;

	int status = cp_evol ( &params, R, t_end, STEP, req_cp, cp_m ) ;
	if ( status != 0 )
		exit (EXIT_FAILURE) ;

	/* free memory for matrices */
 	gsl_matrix_free(cp_m) ;	

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
