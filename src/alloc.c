/*
   -	alloc.c --
   -		This file defines memory allocators.  See alloc (3).
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
   .
   .	$Revision: 1.32 $ $Date: 2014/01/24 22:03:37 $
 */

#include "unix_defs.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "alloc.h"

static int init;
static void alloc_init(void);
static void clean(void);

/* This counter records the number of times an allocator
 * has been called.  It helps to sort output for debugging. */
static unsigned c;

/* Where to send diagnostic output */
static FILE *diag_out;

/* File and line at which to induce pretend memory failure */
static char *fail_fnm;
static int fail_line;

/* Initialize this interface, when process starts */
static void alloc_init(void)
{
    char *s;
    int od;

    if (init) {
	return;
    }
    s = getenv("MEM_DEBUG");
    if (s) {
	if (sscanf(s, "%d", &od) == 1) {
	    diag_out = fdopen(od, "w");
	} else {
	    diag_out = fopen(s, "w");
	}
	if ( !diag_out ) {
	    perror("MEM_DEBUG set but unable to open diagnostic memory file");
	}
	atexit(clean);
    }
    s = getenv("MEM_FAIL");
    if (s) {
	fail_fnm = malloc(strlen(s) + 1);
	if (sscanf(s, "%[^:]:%d", fail_fnm, &fail_line) != 2) {
	    fprintf(stderr, "Could not get failure spec from %s\n", s);
	    free(fail_fnm);
	}
    }
    init = 1;
}

/* Clean up when process exits */
void clean()
{
    if (diag_out) {
	fclose(diag_out);
    }
    if (fail_fnm) {
	free(fail_fnm);
    }
}

/* See alloc (3) */
void *Tkx_Malloc(size_t sz, char *fnm, int ln)
{
    void *m;

    if ( !init ) {
	alloc_init();
    }
    m = malloc(sz);
    if (fail_fnm && (ln == fail_line) && strcmp(fail_fnm, fnm) == 0) {
	return NULL;
    }
    if (m && diag_out) {
	fprintf(diag_out, "%p (%09x) allocated at %s:%d\n", m, ++c, fnm, ln);
    }
    return m;
}

/* See alloc (3) */
void *Tkx_Calloc(size_t n, size_t sz, char *fnm, int ln)
{
    void *m;

    if ( !init ) {
	alloc_init();
    }
    if (fail_fnm && (ln == fail_line) && strcmp(fail_fnm, fnm) == 0) {
	return NULL;
    }
    m = calloc(n, sz);
    if (m && diag_out) {
	fprintf(diag_out, "%p (%09x) allocated at %s:%d\n", m, ++c, fnm, ln);
    }
    return m;
}

/* See alloc (3) */
void *Tkx_ReAlloc(void *m, size_t sz, char *fnm, int ln)
{
    void *m2;

    if ( !init ) {
	alloc_init();
    }
    if (fail_fnm && (ln == fail_line) && strcmp(fail_fnm, fnm) == 0) {
	return NULL;
    }
    m2 = realloc(m, sz);
    if (m2 && diag_out) {
	if (m2 != m) {
	    if (m) {
		fprintf(diag_out, "%p (%09x) freed by realloc at %s:%d\n",
			m, ++c, fnm, ln);
	    }
	    fprintf(diag_out, "%p (%09x) allocated by realloc at %s:%d\n",
		    m2, ++c, fnm, ln);
	} else {
	    fprintf(diag_out, "%p (%09x) reallocated at %s:%d\n",
		    m, ++c, fnm, ln);
	}
    }
    return m2;
}

/* See alloc (3) */
void Tkx_Free(void *m, char *fnm, int ln)
{
    if ( !init ) {
	alloc_init();
    }
    if (diag_out) {
	fprintf(diag_out, "%p (%09x) freed at %s:%d\n", m, ++c, fnm, ln);
    }
    free(m);
}

