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

/* funcs.h */

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define WS_SZ	1000  			/* size of the integration workspace */

#define POW_2	gsl_pow_2  		/* square */
#define VECTOR	gsl_vector_get		/* vector_get */


struct f_params { double omega_c ; double beta ; double Omega ;
	double omega_1; double alpha; } ;
int assign_p ( void*, double*, double*, double*, double* ) ;

struct ent_pars { const gsl_vector* rhoeq ; gsl_matrix* red_m ; } ;

/* Re_gsc functions */
double k_func ( double, void* ) ;
double k_func_1 ( double , void* ) ;
double k_func_2 ( double , void* ) ;
double k_func_3 ( double , void* ) ;
double k_func_4 ( double , void* ) ;
int first_int ( double*, double*, void* ) ;
int second_int ( double* , double* , void* ) ;
int third_int ( double* , double* , void* ) ;
int fourth_int ( double* , double* , void* ) ;
int re_gsc ( double* result, double* error, void* params ) ;
int re_gcs ( double* result, double* error, void* params ) ;

int re_gcc ( void* params, double* val ) ;
int im_gcc ( void* params, double* val, double* error ) ;
int re_gss ( void* params, double* val ) ;
int im_gss ( void* params, double* val, double* error ) ;
int im_gcs ( void* params, double* val ) ;
int im_gsc ( void* params, double* val ) ;
int re_gc0 ( void* params, double* val ) ;
int im_gc0 ( void* params, double* val, double* error ) ;
int re_gs0 ( void* params, double* val, double* error ) ;
int im_gs0 ( void* params, double* val ) ;

/* Img_c0, Img_cc, Img_ss, Img_cs functions */
int expi ( double , double*, double* ) ;
int expi_plus ( double, double*, double* ) ;
double fu ( double, void* ) ;
double ex ( double, void* ) ;

/* Reg_s0 functions */
double fu_inv ( double, void* ) ;
double fu_cau ( double, void* ) ;
double fu_coth ( double, void* ) ;

/* Functions to create matrices */
int red_mat ( gsl_matrix* , void* ) ;
int cp_mat ( gsl_matrix*, void* ) ;
int ham_gen ( gsl_matrix*, const double* ) ;
int mat_write ( gsl_matrix*, char* ) ;
int mat_read ( gsl_matrix*, char* ) ;
int integration ( void*, double* ) ;

/* Functions for evolution */
int generator ( double, const double*, double*, void* ) ;
int jac ( double, const double*, double*, double*, void* ) ;
int evol ( double, gsl_vector*, double, gsl_odeiv2_evolve*,
		gsl_odeiv2_control*, gsl_odeiv2_step*, gsl_odeiv2_system* ) ;
int red_evol ( void*, const double*, double, double, const gsl_vector*,	gsl_matrix* ) ;
int cp_evol ( void*, const double*, double, double, const gsl_vector*,	gsl_matrix* ) ;

/* Stationarity */
int stationary ( const gsl_matrix* , gsl_vector*  ) ;
double J ( double, double ) ;
double polarization ( void*, double ) ;

/* Entropy */
double entropy_of_state ( const gsl_vector* ) ;
double entropy_production ( const gsl_vector*, const gsl_vector*, const gsl_matrix* ) ;
double entropy_prod_time  ( double, void* ) ;
double entropy_integ ( double, const gsl_vector*, gsl_matrix* ) ;

/* DC current */
double current ( const gsl_vector* ) ;
double current_temp ( double ) ;

/* total current */
double tot_current ( const gsl_vector*, void* ) ;

/* writing operations */
int save_integrals ( void* ) ;
int save_matrices ( void* ) ;
int write_red_curr_T ( void* ) ;
int write_cp_curr_T ( void* ) ;
int write_red_curr_O ( void* ) ;
int write_cp_curr_O ( void* ) ;


/* currents */
double current_red_om ( double, void* ) ;
double current_red_T ( double, void* ) ;
double current_cp_om ( double, void* ) ;
double current_cp_T ( double, void* ) ;

/* positivity */
double r0_dot ( const gsl_matrix* , const gsl_vector*  ) ;
int bloch_vector ( gsl_vector* , double , double , double ) ;
int polars ( double* , double* , double* , const gsl_vector*  ) ;
int sample ( const gsl_matrix* , unsigned int ) ;

/* posit */
gsl_matrix_complex* state_matrix ( double, double ) ;
gsl_matrix_complex* rotate ( gsl_matrix_complex*, double ) ;
gsl_matrix_complex* D_psi ( gsl_matrix_complex* ) ;

/* pos_study */
double delta ( double, double*, double, double, double ) ;
int read_integs ( double* ) ;

/* entropy production with different parameters */
int ent_prod_fixed_t ( double, double ) ;
int ent_prod_red_evol ( void* , const double* , double , double , const gsl_vector* , gsl_matrix* , double ) ;
char* dtoascii ( double, int ) ;
