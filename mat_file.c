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
 *       Filename:  mat_file.c
 *
 *    Description:  Write and read the generator matrices on files
 *
 *        Version:  1.0
 *        Created:  03/05/2014 13:52:15
 *       Revision:  none
 *        License:  BSD
 *
 *         Author:  Giuseppe Argentieri (ga), giuseppe.argentieri@ts.infn.it
 *   Organization:  Università degli Studi di Trieste
 *
 * 
 */

#include <gsl/gsl_matrix.h>
#include "funcs.h"

/* 
 *      FUNCTION  
 *         Name:  mat_write
 *  Description:  Save the generator matrix in a file.
 * 
 */
int mat_write ( gsl_matrix* mat, char* name )
{
	FILE* f = fopen ( name, "w" ) ;
	int status = gsl_matrix_fprintf ( f, mat, "%.6f" ) ;
	fclose (f) ;	

	return status;
}		/* -----  end of function mat_write  ----- */


/* 
 *      FUNCTION  
 *         Name:  mat_read
 *  Description:  Read the Redfield matrix from a file 
 * 
 */
int mat_read ( gsl_matrix* mat, char* name )
{
	FILE* f = fopen ( name, "r" ) ;
	int status = gsl_matrix_fscanf ( f, mat ) ;
	fclose (f) ;
	return status;
}		/* -----  end of function mat_read  ----- */
