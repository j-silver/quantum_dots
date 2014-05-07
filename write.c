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
 *       Filename:  write.c
 *
 *    Description:  Write the data into file .dat
 *
 *        Version:  1.0
 *        Created:  07/05/2014 22:46:50
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "funcs.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>

/* 
 *      FUNCTION  
 *         Name:  save_integrals
 *  Description:  
 * 
 */
int save_integrals ( void* params )
{
	/* perform the integrals */
  	double integrals[12] ;
	void* p = (void*) params ;
	int status = integration ( p, integrals ) ;

	/* write them into a file */
	FILE* f_integ = fopen ( "INTEGRALS.dat", "w+" ) ;
	int i ;
	for ( i = 0 ; i < 12 ; i++ )
		fprintf( f_integ, "integrals[%d]: %.6f\n", i, integrals[i] ) ;
	fclose( f_integ ) ;

	return status;
}		/* -----  end of function save_integrals  ----- */


/* 
 *      FUNCTION  
 *         Name:  save_matrices
 *  Description:  
 * 
 */
int save_matrices ( void* params )
{
	double integrals[12] ;

	FILE* f_integ = fopen ( "INTEGRALS.dat", "r" ) ;
	int i ;
	for ( i = 0 ; i < 12 ; i++ )
		fscanf ( f_integ, "%*s %lf", &integrals[i] ) ;
	fclose ( f_integ ) ;
	
	void* p = (void *) params ;

	/* create the Redfield matrix and save it into a file */
	gsl_matrix* red_matrix = gsl_matrix_calloc ( 4, 4 ) ;
	int status1 = red_mat ( red_matrix, integrals, p ) ;
	int status2 = mat_write ( red_matrix, "REDFIELD_MATRIX" ) ;

	/* create the CP matrix and save it into a file */
	gsl_matrix* cp_matrix = gsl_matrix_calloc ( 4, 4 ) ;
	int status3 = cp_mat ( cp_matrix, integrals, p ) ;
	int status4 = mat_write ( cp_matrix, "CP_MATRIX" ) ;

	return status1+status2+status3+status4 ;
}		/* -----  end of function save_matrices  ----- */


