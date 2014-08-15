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
 *       Filename:  3dposit.c
 *
 *    Description:  3D portion of the Bloch sphere where positivity is violated
 *
 *        Version:  1.0
 *        Created:  30/05/2014 17:15:44
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include "s2plot.h"
#include <stdio.h>


/* 
 *      FUNCTION  
 *         Name:  max
 *  Description:  Pick the max value of the dr0/dt derivative 
 * 
 */
double max()
{
	FILE* f = fopen ( "POS_VIOLATIONS", "r" ) ;
	double max = 0 ; double d ;
	while ( fscanf( f, "%*f %*f %*f %lf\n", &d ) == 1 )
		( max < d ) ? max = d : max ;
	fclose (f) ;

	return max;
}		/* -----  end of function max  ----- */


#include	<stdlib.h>

/* 
 *      FUNCTION  
 *         Name:  main
 *  Description:  
 * 
 */
int main ( int argc, char *argv[] )
{
	XYZ      p ;
	COLOUR col ;
	double   D, M ;
	
	M = max() ;
	
	s2opend("/?", argc, argv) ;
	s2swin(-1,1,-1,1,-1,1) ;
 	s2box("BCDETOQ",0,0,"BCDETOQ",0,0,"BCDETOQ",0,0) ;
 
	s2lab("X", "Y", "Z", "Positivity violations") ;

	/* draw blue axes */
	ns2line(-1,0,0,1,0,0, 0,0,1) ;
	ns2line(0,-1,0,0,1,0, 0,0,1) ;
	ns2line(0,0,-1,0,0,1, 0,0,1) ;

	/* white background */
 	ss2sbc(1,1,1) ;

	/* blue foreground */
 	ss2sfc(0,0,1) ;

	FILE* f = fopen ( "POS_VIOLATIONS", "r" ) ;
	while ( fscanf (f, "%lf %lf %lf %lf", &p.x, &p.y, &p.z, &D) == 4 )
	{
		col.r = D/M ;
		ns2vpoint (p, col) ;
	}

	fclose (f) ;

	s2show(1);

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
