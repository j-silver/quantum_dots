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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>

const double omega_c = 1000 ;                  /* critical ohmic frequency */
const double T = 0.1 ;                           /* temperature */
const double Omega = 2.0 ;                     /* pumping frequency */
const double D = 1.0 ;                         /* pumping amplitude */
const double alpha = 5e-3 ;                    /* coupling strength */

const double t_end = 300 ;                     /* time end */
const double step = 0.01 ;                     /* time step */

const double r[] = { 1, 0, -0.894, -0.447 } ;  /* initial state: |z,-> */
const double om[] = { 0, 0, 0, 1 } ;

const int dynamics = 0 ;                       /* 0 = Redfield 
						  1 = CP 
						  2 = Hamiltonian */

int main ( int argc, char* argv[] )
{
	gsl_ieee_env_setup () ;			  /* read GSL_IEEE_MODE */

	double beta = 1/T ;                       /* Boltzmann factor: beta */
	double omega_1 = sqrt(Omega*Omega+D*D) ;  /* omega' */

	/* set the parameters to be passed to the functions */
	struct f_params params;
	params.omega_c = omega_c ;
	params.beta = beta ;
	params.Omega = Omega ;
	params.omega_1 = omega_1 ;
	params.alpha = alpha ;
 
	/* perform the integrals */
  	double integrals[12] ;
  	int status1 = integration ( &params, integrals ) ;

	/* set the Redfield matrix and save onto a file */
	gsl_matrix* red_matrix = gsl_matrix_calloc ( 4, 4 ) ;
	int status2 = red_mat ( red_matrix, integrals, &params ) ;
	int status3 = mat_write ( red_matrix, "REDFIELD_MATRIX" ) ;

	/* read the Redfield matrix from a file */
	gsl_matrix* red_m = gsl_matrix_calloc ( 4, 4 ) ;
	int status7 = mat_read ( red_m, "REDFIELD_MATRIX" ) ;

	/* set the CP matrix and save onto a file */
	gsl_matrix* cp_matrix = gsl_matrix_calloc ( 4, 4 ) ;
	int status4 = cp_mat ( cp_matrix, integrals, &params ) ;
	int status5 = mat_write ( cp_matrix, "CP_MATRIX" ) ;

	/* read the CP matrix from a file */
	gsl_matrix* cp_m = gsl_matrix_calloc ( 4, 4 ) ;
	int status6 = mat_read ( cp_m, "CP_MATRIX" ) ;
	
	/* Hamiltonian generator */
	gsl_matrix* H = gsl_matrix_calloc ( 4, 4 ) ;
	int status8 = ham_gen ( H, om ) ;

	void* PARS ;
	switch ( dynamics ) {
		case 0: 
			printf("Redfield dynamics\n") ;	
			PARS = (void*) red_m ;
			break;
		case 1:
			printf("CP dynamics\n") ;
			PARS = (void*) cp_m ;
			break;
		case 2:
			printf("Hamiltonian dynamics\n") ;
			PARS = (void*) H ;
			break;
	}

	/* Find the stationary state by linear system */
	gsl_vector* req = gsl_vector_calloc ( 4 ) ;
	gsl_vector_set ( req, 0, 1 ) ;
	int status9 = stationary ( (const gsl_matrix*) PARS , req ) ;
	printf("Stationary state: (%.1f , %.6f , %.6f , %.6f )\n", 
			gsl_vector_get(req,0), gsl_vector_get(req,1),
			gsl_vector_get(req,2), gsl_vector_get(req,3) ) ;

	/* 
	 *
	 * Initializing the system
	 *
	 */
	gsl_odeiv2_system sys = { generator, jac, 4, PARS } ;
	gsl_odeiv2_step* s = NULL ; 
	gsl_odeiv2_control* c = NULL ;
	gsl_odeiv2_evolve* e = NULL ;

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 ) ;

	/* Setting the step control: abserr=1e-6, relerr=1e-3 */
	c = gsl_odeiv2_control_standard_new ( 1e-6, 1e-3, 1, 1 ) ;

	/* Allocating the space for evolution function */
	e = gsl_odeiv2_evolve_alloc ( 4 ) ;

	/* 
	 *
	 *  evolving the system 
	 * 
	 */
	double t = 0 ; 
	/* setting the initial vector */
	gsl_vector* init = gsl_vector_calloc(4) ;
	gsl_vector_set ( init, 0, r[0] ) ;
	gsl_vector_set ( init, 1, r[1] ) ;
	gsl_vector_set ( init, 2, r[2] ) ;
	gsl_vector_set ( init, 3, r[3] ) ;

	/* opening the files */
	FILE* f = fopen ( "EVOLUTION.dat", "w" ) ;
	FILE* g = fopen ( "ENTROPY.dat", "w" ) ;

	while ( t < t_end )
	{
		evol ( t, init, step, e, c, s, &sys ) ;
		fprintf ( f, "%.3f %.6f %.6f %.6f\n", t, gsl_vector_get(init,1),
				gsl_vector_get(init,2), gsl_vector_get(init,3) ) ;
		fprintf ( g, "%.3f %.6f\n", t, entropy_production( init , req, 
					(const gsl_matrix*) PARS ) ) ;
		t += step ;
	}

	/*  close the files */
	fclose (f) ;
	fclose (g) ;

	/* free memory for evolution */
	gsl_odeiv2_evolve_free (e) ;
	gsl_odeiv2_control_free (c) ;
	gsl_odeiv2_step_free (s) ;

	/* polarization */
	double D30 = gsl_matrix_get(red_m,3,0) ;
	double D33 = gsl_matrix_get(red_m,3,3) ;
	double pol = -D30/D33 ;
	printf("n-Polarization -D30/D33: %.6f\n", pol ) ;

	return status1 + status2 + status3 + status4 + status5 + status6 + status7 
		+ status8 + status9 ;
}

