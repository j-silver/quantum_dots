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

/* main.c */

#include "funcs.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>

const double omega_c = 1000 ;                   /* critical ohmic frequency */
const double T = 0.1 ;                          /* temperature */
const double D = 1 ;                            /* pumping amplitude (GHz) */
const double Omega = 2 ;                        /* pumping frequency (GHz) */
const double alpha = 5e-3 ;                     /* coupling strength */

const double gamma0 = 0.05 ;                    /* energy hopping between sites
						 * in MeV                       */

const double t_end = 300 ;                      /* time end */
const double step = 0.01 ;                      /* time step */

const double r[] = { 1, 0, -0.894, -0.447 } ;   /* initial state: |z,-> */


int main ( int argc, char* argv[] )
{
	gsl_ieee_env_setup () ;			/* read GSL_IEEE_MODE */

	double beta = 1/T ;                     /* Boltzmann factor: beta */
	double omega_1 = gsl_hypot(Omega,D) ;   /* omega' */


	struct f_params params;
	params.omega_c = omega_c ;
	params.beta = beta ;
	params.Omega = Omega ;
	params.omega_1 = omega_1 ;
	params.alpha = alpha ;
 
	int i ;                                 /* counter for the for loops */

	int status1 = save_integrals ( &params ) ;

	int status2 = save_matrices ( &params ) ;
	
	/* read the Redfield matrix from a file */
	gsl_matrix* red_m = gsl_matrix_calloc ( 4, 4 ) ;
	int status3 = mat_read ( red_m, "REDFIELD_MATRIX" ) ;

	/* read the CP matrix from a file */
	gsl_matrix* cp_m = gsl_matrix_calloc ( 4, 4 ) ;
	int status4 = mat_read ( cp_m, "CP_MATRIX" ) ;
	
	/* Hamiltonian generator */
	const double om[] = { 0 , 0 , 0 , omega_1/2 } ;
	gsl_matrix* H = gsl_matrix_calloc ( 4, 4 ) ;
	int status5 = ham_gen ( H, om ) ;


	/* 
	 * Find the stationary state by solving the linear systems
	 * and the associated stationary currents
	 */
	gsl_vector* req_red = gsl_vector_calloc ( 4 ) ;
	gsl_vector* req_cp  = gsl_vector_calloc ( 4 ) ;
	gsl_vector_set ( req_red, 0, 1 ) ;
	gsl_vector_set ( req_cp, 0, 1 ) ;

	printf("REDFIELD DYNAMICS\n") ;
	int status6 = stationary ( (const gsl_matrix*) (void*) red_m , req_red ) ;
	printf("Stationary state: ( %.1f , %.6f , %.6f , %.6f )\n", 
			gsl_vector_get(req_red,0), gsl_vector_get(req_red,1),
			gsl_vector_get(req_red,2), gsl_vector_get(req_red,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.6f\n\n",
			-gsl_vector_get(req_red,3)*Omega/omega_1) ;

	printf("CP DYNAMICS\n") ;
	int status7 = stationary ( (const gsl_matrix*) (void*) cp_m , req_cp ) ;
	printf("Stationary state: ( %.1f , %.6f , %.6f , %.6f )\n", 
			gsl_vector_get(req_cp,0), gsl_vector_get(req_cp,1),
			gsl_vector_get(req_cp,2), gsl_vector_get(req_cp,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.6f\n\n",
			-gsl_vector_get(req_cp,3)*Omega/omega_1) ;

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

	/* Setting the step control: abserr=1e-6, relerr=1e-3 */
	gsl_odeiv2_control* r_c = gsl_odeiv2_control_standard_new ( 1e-6, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* r_e = gsl_odeiv2_evolve_alloc ( 4 ) ;


	/* 
	 *
	 * Initializing the system for CP dynamics
	 *
	 */
	gsl_odeiv2_system cp_sys = { generator, jac, 4, (void*) cp_m } ;

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	gsl_odeiv2_step* cp_s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 ) ;

	/* Setting the step control: abserr=1e-6, relerr=1e-3 */
	gsl_odeiv2_control* cp_c = gsl_odeiv2_control_standard_new ( 1e-6, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* cp_e = gsl_odeiv2_evolve_alloc ( 4 ) ;


	/* 
	 *
	 *  evolving the systems
	 * 
	 */
	double t = 0 ; 

	/* setting the initial vector */
	gsl_vector* init_red = gsl_vector_calloc(4) ;
	gsl_vector* init_cp = gsl_vector_calloc(4) ;
	for ( i = 0 ; i < 3 ; i++ )
	{
		gsl_vector_set ( init_red, i, r[i] ) ;
		gsl_vector_set ( init_cp, i, r[i] ) ;
	}

	/* opening the files */
	FILE* f_red = fopen ( "RED-EVOLUTION.dat", "w" ) ;
	FILE* f_cp = fopen ( "CP-EVOLUTION.dat", "w" ) ;
	FILE* g_red = fopen ( "RED-ENTROPY.dat", "w" ) ;
	FILE* g_cp = fopen ( "CP-ENTROPY.dat", "w" ) ;
	FILE* h_red = fopen ( "RED-CURRENT.dat", "w" ) ;
	FILE* h_cp = fopen ( "CP-CURRENT.dat", "w" ) ;

	/* writing column heads */
	fprintf ( f_red, "t a(t) b(t) c(t)\n" ) ;
	fprintf ( f_cp,  "t a(t) b(t) c(t)\n" ) ;
	fprintf ( g_red, "t entropy\n" ) ;
	fprintf ( g_cp,  "t entropy\n" ) ;
	fprintf ( h_red, "t I/I0\n" ) ;
	fprintf ( h_cp,  "t I/I0\n" ) ;

	/* writing data */
	while ( t < t_end )
	{
		evol ( t, init_red, step, r_e, r_c, r_s, &red_sys ) ;
		fprintf ( f_red, "%.3f %.6f %.6f %.6f\n", t, gsl_vector_get(init_red,1),
				gsl_vector_get(init_red,2), gsl_vector_get(init_red,3) ) ;
		fprintf ( g_red, "%.3f %.6f\n",
				t, entropy_production( init_red, req_red, red_m )) ;
		fprintf ( h_red, "%.3f %.6f\n",
				t, -gsl_vector_get(init_red,3)*Omega/omega_1 ) ;

		evol ( t, init_cp, step, cp_e, cp_c, cp_s, &cp_sys ) ;
		fprintf ( f_cp, "%.3f %.6f %.6f %.6f\n", t, gsl_vector_get(init_cp,1),
				gsl_vector_get(init_cp,2), gsl_vector_get(init_cp,3) ) ;
		fprintf ( g_cp, "%.3f %.6f\n",
				t, entropy_production( init_cp , req_cp, cp_m )) ;
		fprintf ( h_cp, "%.3f %.6f\n",
				t, -gsl_vector_get(init_cp,3)*Omega/omega_1 ) ;

		t += step ;
	}

	/*  close the files */
	fclose (f_red) ; fclose (f_cp) ;
	fclose (g_red) ; fclose (g_cp) ;
	fclose (h_red) ; fclose (h_cp) ;

	/* free memory for evolution */
	gsl_odeiv2_evolve_free (r_e) ; gsl_odeiv2_evolve_free (cp_e) ;
	gsl_odeiv2_control_free (r_c) ;	gsl_odeiv2_control_free (cp_c) ;
	gsl_odeiv2_step_free (r_s) ; gsl_odeiv2_step_free (cp_s) ;

	/* polarization */
	double D30 = gsl_matrix_get(red_m,3,0) ;
	double D33 = gsl_matrix_get(red_m,3,3) ;
	double pol = -D30/D33 ;
	printf("n-Polarization -D30/D33: %.6f\n", pol ) ;

	return status1 + status2 + status3 + status4 + status5 + status6 + status7 ;
}

