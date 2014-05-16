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
#include "initial.h"

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>


/* 
 *      FUNCTION  
 *         Name:  stationary
 *  Description:  Given the dissipator in Bloch form, reduce to a 3x3 problem and store
 *  			the stationary state in the 3x1 vector *X 
 * 			
 * 			M X = 0
 * 
 * 	        	|  0    0    0    0  |  | 1  |    0
 * 		        | M10  M11  M12  M13 |  | X1 |    0
 * 		        | M20  M21  M22  M23 |  | X2 | =  0
 * 			| M30  M31  M32  M33 |  | X3 |    0
 *
 *
 * 			A x = b
 *
 * 			| M11  M12  M13 |  | X1 |   | -M10 |
 * 			| M21  M22  M23 |  | X2 | = | -M20 |
 * 			| M31  M32  M33 |  | X3 |   | -M30 |
 */
int stationary ( const gsl_matrix* M, gsl_vector* stat_state )
{
	/* Store space for the stationary state */
	gsl_vector* req = gsl_vector_calloc ( 4 ) ;
	gsl_vector_set ( req, 0, 1 ) ;

	/* Copy the dissipator matrix in a temporary local matrix m
	 * (because the algorithm destroys it...) */
	gsl_matrix* m = gsl_matrix_calloc ( 4, 4 ) ;
	gsl_matrix_memcpy ( m, M ) ;

	/* Create a view of the spatial part of vector req */
	gsl_vector_view x = gsl_vector_subvector ( req, 1, 3 ) ;

	/* Create a submatrix view of the spatial part of m and a vector view
	 * of the spatial part of the 0-th column, which goes into -b in the system
	 * A X = b */
	gsl_matrix_view A = gsl_matrix_submatrix ( m, 1, 1, 3, 3 ) ;
	gsl_vector_view b = gsl_matrix_subcolumn ( m, 0, 1, 3 ) ;
	int status1 = gsl_vector_scale ( &b.vector, -1.0 ) ;	

	/* Solve the system A x = b using Householder transformations.
	 * Changing the view x of req => also req is changed, in the spatial part */
	int status2 = gsl_linalg_HH_solve ( &A.matrix, &b.vector, &x.vector ) ;

	/* Set the returning value for the state stat_state */
	*stat_state = *req ;

	/* Free memory */
	gsl_matrix_free(m) ;
	
	return status1 + status2 ;
}		/* -----  end of function stationary  ----- */




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
 
	unsigned int i ;                        /* counter for the for loops */

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

