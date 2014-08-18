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
 *       Filename:  dtoascii.c
 *
 *    Description:  convert a double in a char array
 *
 *        Version:  1.0
 *        Created:  16/08/2014 16:25:00
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), gius.argentieri@gmail.com
 *   Organization:  
 *
 * 
 */

#include <stdlib.h>
#include "funcs.h"

char* fcvt ( double, int, int*, int* );

/* 
 *      FUNCTION  
 *         Name:  dtoascii
 *  Description:  convert a double x to a char array with ndigit 
 *  			digits after the decimal point
 * 
 */
char* dtoascii ( double x, int ndigit )
{
	int decpt, sgn;                         /* decimal point pos. and sign */
	char* s = fcvt(x, ndigit, &decpt, &sgn); /* string w/o d.p. or sign */
	int sz = 0;
	while ( s[sz] != '\0' )
		sz++;                           /* size of string */

	int SZ = sz;                   		/* real size */
	if ( decpt > 0 )                       
		SZ++; 
	if ( decpt <= 0 )
		SZ = SZ - decpt + 2;
	if ( sgn != 0 )                         /* negative sign */
		SZ++;

	char* str = calloc((size_t) SZ, 1);     /* final string */

	int i, first;                     	/* first digit or d.p. */
	if ( sgn != 0 ) {
		if ( decpt <= 0 ) {
			str[0] = '-';
			str[1] = '0';
			str[2] = '.';
			for ( i = 3; i <= (-decpt+2); i++ )
				str[i] = '0';
			first = -decpt+3;
		}
		else {
			str[0] = '-';
			first = 1;
			str[first+decpt] = '.';
		}
	}
	else if ( decpt <= 0 ) {
			str[0] = '0';
			str[1] = '.';
			for ( i = 2; i <= (-decpt+1); i++ )
				str[i] = '0';
			first = -decpt+2;
	}
	else {
		first = 0;
		str[decpt] = '.';
	}
	
	int k;
	if ( sgn == 0 ) {
		for ( k = 0; k < decpt; k++ ) 
			str[first+k] = s[k];
		if ( decpt > 0 )
			for ( k = 1; k <= sz; k++ )
				str[decpt+k] = s[decpt+k-1];
		if ( decpt < 0 )
			for ( k = 0; k <= sz-1; k++ )
				str[first+k] = s[k];
		if ( decpt == 0 )
			for ( k = 1; k <= sz; k++ )
				str[1+k] = s[k-1];
	}
	else {
		for ( k = 0; k < decpt; k++ )
			str[first+k] = s[k];
		if ( decpt > 0 )
			for ( k = 1; k <= sz; k++ )
				str[1+decpt+k] = s[decpt+k-1];
		if ( decpt < 0 )
			for ( k = 0; k <= sz-1; k++ )
				str[first+k] = s[k];
		if ( decpt == 0 )
			for ( k = 1; k <= sz; k++ )
				str[2+k] = s[k-1];
	}


	return (str);
}		/* -----  end of function dtoascii  ----- */
