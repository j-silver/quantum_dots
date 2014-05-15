/*
 *
 *
 *       Filename:  test.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  14/05/2014 19:24:54
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_matrix.h>


#include	<stdlib.h>
#include 	"funcs.h"

/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	/* set a unit Bloch vector (pure state) */
	gsl_vector* R = gsl_vector_calloc (4) ;
	gsl_vector_set(R,0,1) ; gsl_vector_set(R,2,-0.894) ; gsl_vector_set(R,3,-0.447) ;

	/* read the Redfield matrix from a file */
	gsl_matrix* red_m = gsl_matrix_calloc ( 4, 4 ) ;
	int status3 = mat_read ( red_m, "REDFIELD_MATRIX" ) ;

	double S = r0_dot ( red_m, R ) ;

	printf("d r0 / dt = %g\n", S ) ;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
