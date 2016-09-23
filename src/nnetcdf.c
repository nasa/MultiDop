/*
   -	nnetcdf.c --
   -		This file defines functions that read
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

#include <string.h>
#include <stdio.h>
#include "alloc.h"
#include "nnetcdf.h"

/* Open a NetCDF file. See nnetcdf (3). */
int NNC_Open(const char *file_nm, jmp_buf error_env)
{
    int status;
    int ncid;

    if ((status = nc_open(file_nm, 0, &ncid)) != 0) {
	fprintf(stderr, "Could not open %s. NetCDF error message is: %s\n",
		file_nm ? file_nm : "(NULL)", nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return ncid;
}

/* Return the size of a NetCDF dimension.  See nnetcdf (3). */
size_t NNC_Inq_Dim(int ncid, const char *name, jmp_buf error_env)
{
    int dimid;
    int status;
    size_t len;

    if ((status = nc_inq_dimid(ncid, name, &dimid)) != 0) {
	fprintf(stderr, "Could not find dimension named %s. "
		"NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_dim(ncid, dimid, NULL, &len)) != 0) {
	fprintf(stderr, "Could not retrieve size of %s dimension.  "
		"NetCDF error message is: %s\n", name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return len;
}

/* Return a string from a NetCDF file. See nnetcdf (3). */
char * NNC_Get_String(int ncid, const char *name, jmp_buf error_env)
{
    char *val;		/* Return value */
    int varid;		/* Variable identifier */
    int dimid;		/* Dimension of variable */
    size_t len;		/* Dimension size */
    char *c, *ce;	/* Loop parameters */
    int status;		/* NetCDF function return value */

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_vardimid(ncid, varid, &dimid)) != 0) {
	fprintf(stderr, "Could not get dimension for %s. "
		"NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_dimlen(ncid, dimid, &len)) != 0) {
	fprintf(stderr, "Could not get dimension length for %s. "
		"NetCDF error message is: %s\n", name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !(val = MALLOC(len + 1)) ) {
	fprintf(stderr, "Allocation failed for %s.\n", name);
	longjmp(error_env, NNCDF_ERROR);
    }
    for (c = val, ce = c + len; c < ce; c++) {
	*c = ' ';
    }
    if ((status = nc_get_var_text(ncid, varid, val)) != 0) {
	fprintf(stderr, "Could not get value for %s. NetCDF error message "
		"is: %s\n", name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    *(val + len) = '\0';

    return val;
}

/* Retrieve a character variable from a NetCDF file. See nnetcdf (3). */
char *NNC_Get_Var_Text(int ncid, const char *name, char *cPtr,
	jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !cPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate %d dimension array for %s\n",
		    ndims, name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n", name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s. "
			"NetCDF error message is: %s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(cPtr = MALLOC(sz)) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_text(ncid, varid, cPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s. "
		"NetCDF error message is: %s\n", name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return cPtr;
}

/* Retrieve an unsigned char variable from a NetCDF file. See nnetcdf (3). */
unsigned char * NNC_Get_Var_UChar(int ncid, const char *name,
	unsigned char *uPtr, jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !uPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n", name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s. "
			"NetCDF error message is:%s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(uPtr = MALLOC(sz * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_uchar(ncid, varid, uPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s.  "
		"NetCDF error message is: %s\n", name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return uPtr;
}

/* Retrieve an integer variable from a NetCDF file. See nnetcdf (3). */
int * NNC_Get_Var_Int(int ncid, const char *name, int *iPtr, jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !iPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n", name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s. "
			"NetCDF error message is:%s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(iPtr = MALLOC(sz * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_int(ncid, varid, iPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s.  "
		"NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return iPtr;
}

/* Retrieve an unsigned integer variable from a NetCDF file. See nnetcdf (3). */
unsigned * NNC_Get_Var_UInt(int ncid, const char *name, unsigned *iPtr,
	jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !iPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n", name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s. "
			"NetCDF error message is:%s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(iPtr = MALLOC(sz * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_uint(ncid, varid, iPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s.  "
		"NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return iPtr;
}

/* Retrieve a float variable from a NetCDF file. See nnetcdf (3). */
float * NNC_Get_Var_Float(int ncid, const char *name, float *fPtr,
	jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !fPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s. "
			"NetCDF error message is:%s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(fPtr = MALLOC(sz * sizeof(float))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_float(ncid, varid, fPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s.  "
		"NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return fPtr;
}

/* Retrieve a double variable from a NetCDF file. See nnetcdf (3). */
double * NNC_Get_Var_Double(int ncid, const char *name, double *dPtr,
	jmp_buf error_env)
{
    int varid;		/* Variable identifier */
    int status;

    if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !dPtr ) {
	int ndims;
	static int *dimidPtr;
	int *d, *de;
	size_t sz;

	if ((status = nc_inq_varndims(ncid, varid, &ndims)) != 0) {
	    fprintf(stderr, "Could not get dimension count for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ( !(dimidPtr = REALLOC(dimidPtr, ndims * sizeof(int))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
	if ((status = nc_inq_vardimid(ncid, varid, dimidPtr)) != 0) {
	    fprintf(stderr, "Could not get dimensions for %s. "
		    "NetCDF error message is: %s\n",
		    name, nc_strerror(status));
	    longjmp(error_env, NNCDF_ERROR);
	}
	for (sz = 1, d = dimidPtr, de = d + ndims; d < de; d++) {
	    size_t l;

	    if ((status = nc_inq_dimlen(ncid, *d, &l)) != 0) {
		fprintf(stderr, "Could not get dimension size for %s."
			" NetCDF error message is:%s\n",
			name, nc_strerror(status));
		longjmp(error_env, NNCDF_ERROR);
	    }
	    sz *= l;
	}
	if ( !(dPtr = MALLOC(sz * sizeof(double))) ) {
	    fprintf(stderr, "Could not allocate dimension array for %s\n",
		    name);
	    longjmp(error_env, NNCDF_ERROR);
	}
    }
    if ((status = nc_get_var_double(ncid, varid, dPtr)) != 0) {
	fprintf(stderr, "Could not get value for %s. "
		" NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return dPtr;
}

/* Get a string attribute associated with a NetCDF variable. See nnetcdf (3). */
char * NNC_Get_Att_String(int ncid, const char *name, const char *att,
	jmp_buf error_env)
{
    int varid;
    int status;
    size_t len;
    char *val;

    if (strcmp(name, "NC_GLOBAL") == 0) {
	varid = NC_GLOBAL;
    } else if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_attlen(ncid, varid, att, &len)) != 0) {
	fprintf(stderr, "Could not get string length for %s of %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !(val = MALLOC(len + 1)) ) {
	fprintf(stderr, "Allocation failed for %s\n",name);
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_get_att_text(ncid, varid, att, val)) != 0) {
	fprintf(stderr, "Could not get %s for %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    val[len] = '\0';
    return val;
}

/* Get integer attribute associated with a NetCDF variable. See nnetcdf (3). */
int *NNC_Get_Att_Int(int ncid, const char *name, const char *att,
	jmp_buf error_env)
{
    int varid;
    int status;
    size_t len;
    int *i;

    if (strcmp(name, "NC_GLOBAL") == 0) {
	varid = NC_GLOBAL;
    } else if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_attlen(ncid, varid, att, &len)) != 0) {
	fprintf(stderr, "Could not get string length for %s of %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !(i = CALLOC(len, sizeof(int))) ) {
	fprintf(stderr, "Allocation failed for %s\n",name);
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_get_att_int(ncid, varid, att, i)) != 0) {
	fprintf(stderr, "Could not get %s attribute for %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return i;
}

/*
   Get unsigned integer attribute associated with a NetCDF variable.
   See nnetcdf (3).
 */

unsigned *NNC_Get_Att_UInt(int ncid, const char *name, const char *att,
	jmp_buf error_env)
{
    int varid;
    int status;
    size_t len;
    unsigned *i;

    if (strcmp(name, "NC_GLOBAL") == 0) {
	varid = NC_GLOBAL;
    } else if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_attlen(ncid, varid, att, &len)) != 0) {
	fprintf(stderr, "Could not get attribute length for %s of %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !(i = CALLOC(len, sizeof(unsigned))) ) {
	fprintf(stderr, "Allocation failed for %s\n",name);
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_get_att_uint(ncid, varid, att, i)) != 0) {
	fprintf(stderr, "Could not get %s attribute for %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return i;
}

/* Get a float attribute associated with a NetCDF variable. See nnetcdf (3). */
float *NNC_Get_Att_Float(int ncid, const char *name, const char *att,
	jmp_buf error_env)
{
    int varid;
    int status;
    size_t len;
    float *v;

    if (strcmp(name, "NC_GLOBAL") == 0) {
	varid = NC_GLOBAL;
    } else if ((status = nc_inq_varid(ncid, name, &varid)) != 0) {
	fprintf(stderr, "No variable named %s. NetCDF error message is: %s\n",
		name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_inq_attlen(ncid, varid, att, &len)) != 0) {
	fprintf(stderr, "Could not get attribute length for %s of %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    if ( !(v = CALLOC(len, sizeof(float))) ) {
	fprintf(stderr, "Allocation failed for %s\n",name);
	longjmp(error_env, NNCDF_ERROR);
    }
    if ((status = nc_get_att_float(ncid, varid, att, v)) != 0) {
	fprintf(stderr, "Could not get %s attribute for %s."
		" NetCDF error message is: %s\n",
		att, name, nc_strerror(status));
	longjmp(error_env, NNCDF_ERROR);
    }
    return v;
}

