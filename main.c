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

	/* Reg_cs */
	double regcs, regcserr ;
	w_integration ( &regcs, &regcserr, &params ) ;
	printf("Re_gcs: %.6f +- %.6f\n", regcs, regcserr ) ;

	/* Reg_cc */
	double regcc ;
	re_gcc ( &params, &regcc ) ;
	printf("Re_gcc: %.6f\n", regcc) ;

	/* Reg_ss */
	double regss ;
	re_gss ( &params, &regss ) ;	
	printf("Re_gss: %.6f\n", regss) ;

	/* Img_cs */
	double imgcs ;
	im_gcs ( &params, &imgcs ) ;
	printf("Im_gcs: %.6f\n", imgcs) ;

	/* Img_sc */
	double imgsc ;
	im_gsc ( &params, &imgsc ) ;
	printf("Im_gsc: %.6f\n", imgsc) ;

	/* Reg_c0 */
	double regc0 ;
	re_gc0 ( &params, &regc0 ) ;
	printf("Re_gc0: %.6f\n", regc0) ;

	/* Img_c0 */
	double imgc0, imgc0_error ;
	im_gc0 ( &params, &imgc0, &imgc0_error ) ;
	printf("Im_gc0: %.6f +- %.6f\n", imgc0, imgc0_error) ;

	/* Img_ss */
	double imgss, imgss_error ;
	im_gss ( &params, &imgss, &imgss_error ) ;
	printf("Im_gss: %.6f +- %.6f\n", imgss, imgss_error) ;

	/* Img_cc */
	double imgcc, imgcc_error ;
	im_gcc ( &params, &imgcc, &imgcc_error ) ;
	printf("Im_gcc: %.6f +- %.8f\n", imgcc, imgcc_error) ;

	return 0;
}

