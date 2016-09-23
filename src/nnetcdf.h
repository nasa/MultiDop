/*
   -	nnetcdf.h --
   -		This file declares functions that read
   -		NetCDF files.  See nnetcdf (3).
   -	
   .	Copyright (c) 2011, Gordon D. Carrie. All rights reserved.
   .	
   .	Redistribution and use in source and binary forms, with or without
   .	modification, are permitted provided that the following conditions
   .	are met:
   .	
   .	    * Redistributions of source code must retain the above copyright
   .	    notice, this list of conditions and the following disclaimer.
   .
   .	    * Redistributions in binary form must reproduce the above copyright
   .	    notice, this list of conditions and the following disclaimer in the
   .	    documentation and/or other materials provided with the distribution.
   .	
   .	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   .	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   .	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   .	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   .	HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   .	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
   .	TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   .	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   .	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   .	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   .	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   .
   .	Please send feedback to dev0@trekix.net
 */

#ifndef NNCDF_H_
#define NNCDF_H_

#include <setjmp.h>
#include <netcdf.h>

#define NNCDF_ERROR 1

int NNC_Open(const char *, jmp_buf);
size_t NNC_Inq_Dim(int, const char *, jmp_buf);
char *NNC_Get_Var_Text(int, const char *, char *, jmp_buf);
char *NNC_Get_String(int, const char *, jmp_buf);
unsigned char *NNC_Get_Var_UChar(int, const char *, unsigned char *, jmp_buf);
int *NNC_Get_Var_Int(int, const char *, int *, jmp_buf);
unsigned *NNC_Get_Var_UInt(int, const char *, unsigned *, jmp_buf);
float *NNC_Get_Var_Float(int, const char *, float *, jmp_buf);
double *NNC_Get_Var_Double(int, const char *, double *, jmp_buf);
char *NNC_Get_Att_String(int, const char *, const char *, jmp_buf);
int *NNC_Get_Att_Int(int, const char *, const char *, jmp_buf);
unsigned *NNC_Get_Att_UInt(int, const char *, const char *, jmp_buf);
float *NNC_Get_Att_Float(int, const char *, const char *, jmp_buf);

#endif
