/* Copyright (c) 2014, Giuseppe Argentieri <gius.argentieri@gmail.com>

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
 *       Filename:  cp_violations.c
 *
 *    Description:  Entropy prod. violations at different parameters values
 *    			(T, Delta, OMEGA) 
 *
 *        Version:  1.0
 *        Created:  15/08/2014 14:26:42
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), gius.argentieri@gmail.com
 *   Organization:  
 *
 * 
 */

#include "funcs.h"
#include "initial.h"

#include <stdlib.h>

char* fcvt (double, int, int*, int*);

/* 
 *      FUNCTION  
 *         Name:  ent_prod_fixed_t
 *  Description:  Entropy production evolution at fixed T, varying 
 *  		  the ratio Omega/Delta
 * 
 */
int ent_prod_fixed_t ( double Te, double x )
{
	static int status = 0;

	double De = 1;
	double O = x*De;
	double omega_1 = gsl_hypot(O, De);

	/* set parameters */
	struct f_params params;
	params.omega_c = omega_c;
	params.beta = 1.0/Te;
	params.Omega = O;
	params.omega_1 = omega_1;
	params.alpha = alpha;

	/* create generator and find stationary state */
	gsl_matrix* red_m = gsl_matrix_calloc(4, 4);
	int s1 = red_mat(red_m, &params);
	gsl_vector* red_stat = gsl_vector_calloc(4);
	int s2 = stationary(red_m, red_stat);

	/* evolve the system and write the entr. prod. */
	int s3 = ent_prod_red_evol(&params, R, t_end, STEP, red_stat,
			red_m, x);

	status += s1+s2+s3;
	return (status);
}		/* -----  end of function ent_prod_fixed_t  ----- */


/* 
 *      FUNCTION  
 *         Name:  ent_prod_red_evol
 *  Description:  
 * 
 */
int ent_prod_red_evol ( void* params, const double r[], double time_end, 
		double step, const gsl_vector* req_red, gsl_matrix* red_m, double x )
{
	unsigned int i;                        /* counter for the for loops */

	/* 
	 *
	 *  evolving the systems
	 * 
	 */
	double t = 0; 

	/* setting the initial vector */
	gsl_vector* init_red = gsl_vector_calloc(4);
	for ( i = 0; i < 4; i++ )
		gsl_vector_set ( init_red, i, r[i] );

	/* 
	 *
	 * Initializing the system for Redfield dynamics
	 *
	 */
	gsl_odeiv2_system red_sys = { generator, jac, 4, (void*) red_m };

	/* Choosing the step function type: Runge-Kutta-Fehlberg (4,5) */	
	/* gsl_odeiv2_step* s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkf45 , 4 ) ; */

	/* Choosing the step function type: Runge-Kutta Cash-Karp (4,5) */	
	gsl_odeiv2_step* r_s = gsl_odeiv2_step_alloc ( gsl_odeiv2_step_rkck , 4 );

	/* Setting the step control: abserr=1e-9, relerr=1e-3 */
	gsl_odeiv2_control* r_c = gsl_odeiv2_control_standard_new ( 1e-9, 1e-3, 1, 1 );

	/* Allocating the space for evolution function */
	gsl_odeiv2_evolve* r_e = gsl_odeiv2_evolve_alloc ( 4 );

	/* opening the file */
	int decpt, sgn;
	char* ratio = fcvt(x, 2, &decpt, &sgn); 
	unsigned int sz = 0;
	while ( ratio[sz] != '\0' )
		sz++;
	char Ratio[sz+1];
	int k;
	for ( k = 0; k < decpt; k++ )
		Ratio[k] = ratio[k];
	if ( decpt >= 0 )
		Ratio[decpt] = '.';
	for ( k = decpt+1; k < sz+1; k++ )
		Ratio[k] = ratio[k-1];
	char S[sz+2];
	S[0] = '0';
	for ( k = 1; k < sz+2; k++ )
		S[k] = Ratio[k-1];

	printf("Ratio: %s\n", S);
	FILE* g_red = fopen ( S, "w" );


	/* writing data */
	while ( t < t_end )
	{
		evol ( t, init_red, step, r_e, r_c, r_s, &red_sys );
		fprintf ( g_red, "%.2f %.9f\n",	t,
				entropy_production( init_red, req_red, red_m ));
		t += step;
	}
	
	/*  close the files */
	fclose (g_red);

	/* free memory for evolution */
	gsl_odeiv2_evolve_free (r_e) ;
	gsl_odeiv2_control_free (r_c) ;
	gsl_odeiv2_step_free (r_s) ;

	return (0);
}		/* -----  end of function ent_prod_red_evol  ----- */



/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	unsigned int i;		/*  counter  */
	double x = 0;		/*  O/D ratio  */

	for ( i = 0; i < 9; i++ )
	{
		x += 0.1;
		ent_prod_fixed_t (T, x);
	}

	x = 1;
	for ( i = 0; i < 9; i++ )
	{
		x += 1;
		ent_prod_fixed_t (T, x);
	}

	x = 10;
	for ( i = 0; i < 9; i++ )
	{
		x += 10;
		ent_prod_fixed_t (T, x);
	}


	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
