/* funcs.h */

#include <math.h>

#define WS_SZ 1000  /* size of the integration workspace */

struct f_params { double omega_c ; double beta ; double Omega ;
	double omega_1; double alpha; } ;

int assign_p ( void*, double*, double*, double*, double* ) ;

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

