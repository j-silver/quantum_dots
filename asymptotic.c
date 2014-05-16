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

#include "funcs.h"
#include "initial.h"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>

int main ( int argc, char* argv[] )
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
	gsl_vector* req_red = gsl_vector_calloc (4) ;
	gsl_vector* req_cp = gsl_vector_calloc (4) ;

	printf("REDFIELD DYNAMICS\n") ;
	int status6 = stationary ( red_m , req_red ) ;
	printf("Stationary state: ( %.1f , %.9f , %.9f , %.9f )\n", 
			VECTOR(req_red,0), VECTOR(req_red,1),
			VECTOR(req_red,2), VECTOR(req_red,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.9f\n\n",
			-VECTOR(req_red,3)*OMEGA/omega_1) ;
	
	printf("CP DYNAMICS\n") ;
	int status7 = stationary ( cp_m , req_cp ) ;
	printf("Stationary state: ( %.1f , %.9f , %.9f , %.9f )\n", 
			VECTOR(req_cp,0), VECTOR(req_cp,1),
			VECTOR(req_cp,2), VECTOR(req_cp,3) ) ;
	printf("Stationary normalized (I/I0) DC current: %.9f\n\n",
			-VECTOR(req_cp,3)*OMEGA/omega_1) ;

	/* Save the stationary states into files */
	FILE* f = fopen ( "RED_STATIONARY.dat", "w" ) ;
	gsl_vector_fprintf ( f, req_red, "%.9f" ) ;
	FILE* g = fopen ( "CP_STATIONARY.dat", "w" ) ;
	gsl_vector_fprintf ( g, req_cp, "%.9f" ) ;

	fclose (f) ; fclose (g) ;

	/* polarization */
	double D30 = gsl_matrix_get(cp_m,3,0) ;
	double D33 = gsl_matrix_get(cp_m,3,3) ;
	double pol = -D30/D33 ;
	printf("n-Polarization of the CP dynamics -D30/D33: %.9f\n", pol ) ;

	/* free memory for matrices */
	gsl_matrix_free(red_m) ;
 	gsl_matrix_free(cp_m) ;	

	return status1 + status2 + status3 + status4 + status5 + status6 + status7 ;
}

