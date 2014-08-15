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
 *       Filename:  pos_study.c
 *
 *    Description:  Study of the positivity violations of the Redfield dynamics
 *    			on pure states
 *
 *        Version:  1.0
 *        Created:  20/07/2014 18:16:03
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include	"funcs.h"
#include	"initial.h"
#include	<stdlib.h>
#include	<math.h>

/* 
 *      FUNCTION  
 *         Name:  Delta
 *  Description:  function Delta(alpha) that measures positivity of
 *  			<phi|D[|psi><psi|]|phi>
 * 
 */
double delta ( double a, double integrals[], double O, double d, double o_1 )
{
	/* define integrals */
	double RS0 = integrals[10] ;
	double RCC = integrals[0]  ;
	double RSS = integrals[2]  ;
	double ICS = integrals[7]  ;
	double ISC = integrals[5]  ;
	double RC0 = integrals[8]  ;
	double RSC = integrals[4]  ;
	double RCS = integrals[6]  ;
	double ISS = integrals[3]  ;
	double ICC = integrals[1]  ;

	double den   = (1+a*a)*(1+a*a) ;
	double a_1 = a*a-1 ;
	double ap1 = a*a+1 ;

	double A ;
	A = (1/den)*(RS0*(-2*(O/o_1)*a*a_1) +RCC*(a_1*a_1) +RSS*(-(d/o_1)*a_1*a_1)
		+ICS*a_1*ap1 +ISC*a_1*ap1*d/o_1 +RC0*(-4*a*a*O/o_1) 
		+RSC*(-2*a*ap1) + RCS*(-(d/o_1)*2*a*ap1) +ISS*(2*a*(-a_1)) 
		+ICC*(d/o_1)*(-a_1)*2*a) ;

	return (A);
}		/* -----  end of function delta  ----- */


/* 
 *      FUNCTION  
 *         Name:  read_integs
 *  Description:  Read integrals from disk 
 * 
 */
int read_integs ( double integs[] )
{
	/* read integrals from file INTEGRALS.dat */
	FILE* f_integ = fopen ( "INTEGRALS.dat", "r" ) ;
	int i = 0 ;
	while ( fscanf (f_integ, "%*s %lf\n", &integs[i]) == 1 )
		i++ ;
	fclose ( f_integ ) ;
	
	return (0);
}		/* -----  end of function read_integs  ----- */


/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	double integrals[12] ;
	read_integs( integrals ) ;

	double O = OMEGA ;
	double o_1 = sqrt(O*O+D*D) ;

	/* write Delta(alpha) on file */
	FILE* dalpha = fopen ( "D_ALPHA.dat", "w+" ) ;

	double a ;
	a = 0 ;
	while ( a < 10 )
	{
		fprintf(dalpha, "%.4f %.9f\n", a, delta(a, integrals, O, D, o_1)) ;
		a+=.0001 ;
	}

	fclose ( dalpha ) ;

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
