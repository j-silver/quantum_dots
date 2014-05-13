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
 *       Filename:  current_tdel.c
 *
 *    Description:  Stationary current as a function of T/D
 *
 *        Version:  1.0
 *        Created:  10/05/2014 19:52:35
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_ieee_utils.h>
#include <stdlib.h>
#include <stdio.h>
#include "funcs.h"

const double omega_c = 1000 ;                   /* critical ohmic frequency */
const double D = 1 ;                            /* pumping amplitude (GHz) */
const double Omega = 2 ;                        /* pumping frequency (GHz) */
const double alpha = 5e-3 ;                     /* coupling strength */


/* 
 *      FUNCTION  
 *         Name:  current_red_T
 *  Description:  Given the physical parameters and the ratio x = T/Delta, 
 *  		  determine the generators matrix in Redfield case
 *  		  and the stationary current.
 *
 */
double current_red_T ( double x, void* params )
{
	struct f_params* pars = (struct f_params*) params ;

	/* Calculate beta and set into pars */
	double T = x*D ;
	pars->beta = 1.0/T ;

	/* Calculate the matrix */
	gsl_matrix* m = gsl_matrix_calloc(4,4) ;	
	int s1 = red_mat ( m, pars ) ;

	/* Calculate the stationary current */
	gsl_vector* stat_state = gsl_vector_calloc (4) ;
	int s2 = stationary ( m , stat_state ) ;
	
	if ( (s1 + s2) != 0 )
		exit(EXIT_FAILURE) ;

	double curr = -gsl_vector_get(stat_state,3)*Omega/pars->omega_1 ;

	gsl_matrix_free(m) ;

	return curr ;
}		/* -----  end of function current_red_T  ----- */


/* 
 *      FUNCTION  
 *         Name:  current_cp_T
 *  Description:  
 * 
 */
double current_cp_T ( double x, void* params )
{
	struct f_params* pars = (struct f_params*) params ;

	/* Calculate beta and set into pars */
	double T = x*D ;
	pars->beta = 1.0/T ;	

	/* Calculate the matrix */
	gsl_matrix* m = gsl_matrix_calloc(4,4) ;	
	int s1 = cp_mat ( m, pars ) ;

	/* Calculate the stationary current */
	gsl_vector* stat_state = gsl_vector_calloc (4) ;
	int s2 = stationary ( m , stat_state ) ;
	
	if ( (s1 + s2) != 0 )
		exit(EXIT_FAILURE) ;

	double curr = -gsl_vector_get(stat_state,3)*Omega/pars->omega_1 ;

	gsl_matrix_free(m) ;

	return curr ;
}		/* -----  end of function current_cp_T  ----- */


/* 
 *      FUNCTION  
 *         Name:  write_red_curr_T
 *  Description:  
 * 
 */
int write_red_curr_T ( void* params )
{
	FILE* f = fopen ( "RED-STAT-CURR-T.dat" , "w" ) ;

	double x = 0.00 ;

	int i ;
	for ( i = 0 ; i < 1000 ; i++ )
	{
		x += 0.001 ;
		fprintf ( f, "%.3f %.9f\n", x, current_red_T (x, params) ) ;
	}

	for ( i = 0 ; i < 900 ; i++ )
	{
		x += 0.01 ;
		fprintf ( f, "%.2f %.9f\n", x, current_red_T (x, params) ) ;
	}

	for ( i = 0 ; i < 900 ; i++ )
	{
		x += 0.1 ;
		fprintf ( f, "%.1f %.9f\n", x, current_red_T (x, params) ) ;
	}

	fclose (f) ;

	return 0;
}		/* -----  end of function write_red_curr_T  ----- */

/* 
 *      FUNCTION  
 *         Name:  write_cp_curr_T
 *  Description:  
 * 
 */
int write_cp_curr_T ( void* params )
{
	FILE* f = fopen ( "CP-STAT-CURR-T.dat" , "w" ) ;

	double x = 0.00 ;

	int i ;
	for ( i = 0 ; i < 1000 ; i++ )
	{
		x += 0.001 ;
		fprintf ( f, "%.3f %.9f\n", x, current_cp_T (x, params) ) ;
	}

	for ( i = 0 ; i < 900 ; i++ )
	{
		x += 0.01 ;
		fprintf ( f, "%.2f %.9f\n", x, current_cp_T (x, params) ) ;
	}

	for ( i = 0 ; i < 900 ; i++ )
	{
		x += 0.1 ;
		fprintf ( f, "%.1f %.9f\n", x, current_cp_T (x, params) ) ;
	}

	fclose (f) ;

	return 0;
}		/* -----  end of function write_cp_curr_T  ----- */



/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	double omega_1 = gsl_hypot(Omega,D) ;   /* omega' */

	struct f_params params;
	params.omega_c = omega_c ;
	params.beta = 1 ;                       /* it will be changed after... */
	params.Omega = Omega ;
	params.omega_1 = omega_1 ;
	params.alpha = alpha ;

	int status1 = write_red_curr_T ( &params ) ;
	int status2 = write_cp_curr_T ( &params ) ;

	return status1+status2 ;
}				/* ----------  end of function main  ---------- */
