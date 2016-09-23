/*
   Apply Leise filter to a 2D array.

   Build like so --
   gfortran -o leise_example leise_filt.f leise_example.c
   Run --
   ./leise_example

 */

#include <stdlib.h>
#include <stdio.h>

static float **falloc2(long, long, const char *);

int main()
{
    int num_x = 7;
    int num_y = 6;
    float **fld;
    int num_z = 1;			/* Always 1 for 2D */
    int num_step = 1;			/* Ignore */
    int k, j, i;

    fld = falloc2(num_y, num_x, "fld");

    /* Initialize fld to 0.0. Then put a box in the middle. */
    for (j = 0; j < num_y; j++) {
	for (i = 0; i < num_x; i++) {
	    fld[j][i] = 0.0;
	}
    }
    for (j = 2; j <= 3; j++) {
	for (i = 2; i <= 4; i++) {
	    fld[j][i] = 4.0;
	}
    }

    /* Print initial field */
    printf("Before filter:\n");
    for (j = 0; j < num_y; j++) {
	for (i = 0; i < num_x; i++) {
	    printf(" %6.3f", fld[j][i]);
	}
	printf("\n");
    }
    printf("\n");

    /* The 1D array in FORTRAN is fld[0] in C */
    t5fltr_(fld[0], &num_x, &num_y, &num_z, &num_step);

    /* Print filtered field */
    printf("Filtered:\n");
    for (j = 0; j < num_y; j++) {
	for (i = 0; i < num_x; i++) {
	    printf(" %6.3f", fld[j][i]);
	}
	printf("\n");
    }
    printf("\n");

    return 0;
}

static float **falloc2(long jmax, long imax, const char *var_nm)
{
    float **dat = NULL;
    long j, i;
    size_t jj, ii;

    /* Make sure casting to size_t does not overflow anything.  */
    if ( jmax <= 0 || imax <= 0 ) {
	fprintf(stderr, "Array dimensions must be positive for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    jj = (size_t)jmax;
    ii = (size_t)imax;
    if ((jj * ii) / jj != ii) {
	fprintf(stderr, "Dimensions [%ld][%ld]  too big "
		"for pointer arithmetic when allocating %s.\n",
		jmax, imax, var_nm);
	exit(EXIT_FAILURE);
    }

    dat = (float **)calloc(jj + 2, sizeof(float **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (float *)calloc(jj * ii + 1, sizeof(float));
    if ( !dat[0] ) {
	fprintf(stderr, "Could not allocate 1st dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    for (j = 1; j <= jmax; j++) {
	dat[j] = dat[j - 1] + imax;
    }
    for (i = 0; i < jmax * imax; i++) {
	dat[0][i] = 0.0;
    }
    return dat;
}
