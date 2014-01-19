/* funcs.h */

#include <math.h>

#define WS_SZ 1000  /* size of the integration workspace */

struct f_params { double omega_c ; double beta ; double Omega ;
	double omega_1; double alpha; } ;

int assign_p ( void*, double*, double*, double*, double* ) ;

/* Re_gcs functions */
double k_integrand ( void*, double k, double w ) ;
double k_func ( double, void* ) ;
double k_integration ( double, void* ) ;
double w_integrand ( double, void* ) ;
int    w_integration ( double*, double*, void* ) ;
struct ext_pars { struct f_params params; double w; } ;int re_gcc ( void* params, double* val ) ;


int re_gss ( void* params, double* val ) ;
int im_gcs ( void* params, double* val ) ;
int im_gsc ( void* params, double* val ) ;
int re_gc0 ( void* params, double* val ) ;
int im_gc0 ( void* params, double* val, double* error ) ;
int im_gss ( void* params, double* val, double* error ) ;
int im_gcc ( void* params, double* val, double* error ) ;

/* Img_c0, Img_cc, Img_cs functios */
int expi ( double , double*, double* ) ;
int expi_plus ( double, double*, double* ) ;
double fu ( double, void* ) ;
double ex ( double, void* ) ;
