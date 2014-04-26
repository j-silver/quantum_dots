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

const double omega_c = 100.0 ;
const double beta = 1.0 ;
const double Omega = 1.0 ;
const double omega_1 = 5.0 ;
const double alpha = 0.01 ;

int main ( int argc, char* argv[] )
{
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

	int i ;
	for ( i = 0 ; i<12 ; i++ )
		printf("integrals[%d]: %.6f\n", i, integrals[i] ) ;

	/* set the Redfield matrix and save onto a file */
	gsl_matrix* red_matrix = gsl_matrix_calloc ( 4, 4 ) ;
	int status2 = red_mat ( red_matrix, integrals, &params ) ;
	int status3 = red_mat_write ( red_matrix, "Redfield_matrix" ) ;

	return status1 + status2 + status3;
}

