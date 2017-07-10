/*
   coverage.c --
   	Store and provide DDA coverage coefficients for various combinations
	of radars.

   $Revision: 1.6 $ $Date: 2014/11/07 17:10:18 $
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "DDA.h"
#include "coverage.h"

/* Number of radars in this analysis */ 
static size_t num_radars_max;

/* Unusable values */ 
#define I_BAD -INT_MAX

/*
   These arrays store coverage values for various combinations
   combinations of radars. Each array has CVG_HASH_SZ elements. Radar
   combination denoted by indeces {i1, i2, ..., in} corresponds
   to offset 2^i1 + 2^i2 ... + 2^in in these arrays. In other words,
   each bit in an index corresponds to a radar.

   Use global function DDA_Get_Radar_Idx to get index for a radar given
   its name.
 */

#define CVG_HASH_SZ 256
static int cvg_opt_bg[CVG_HASH_SZ];
static int cvg_sub_bg[CVG_HASH_SZ];
static int cvg_opt_fil[CVG_HASH_SZ];
static int cvg_sub_fil[CVG_HASH_SZ];
static int cvg_bg[CVG_HASH_SZ];
static int cvg_fil[CVG_HASH_SZ];
static double sseq_trip[CVG_HASH_SZ];

/* Local functions */
static void set_cvg_opt_bg(int, int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_cvg_sub_bg(int, int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_cvg_opt_fil(int, int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_cvg_sub_fil(int, int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_cvg_bg(int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_cvg_fil(int, char [MAX_RADARS][PATH_LEN_MAX]);
static void set_sseq_trip(double, char [MAX_RADARS][PATH_LEN_MAX]);
static int cvg_idx_i(int, int *);
static unsigned int cvg_idx_s(int, char [MAX_RADARS][PATH_LEN_MAX]);

/*
   This function reads coverage parameters from the file named fl_nm.
   See sample calc_params.dda file for explanation of input.
   Local functions are set.  Return value is 1/0 for success/failure.
 */

int Coverage_Read(const char *fl_nm, size_t num_radars_in)
{
    FILE *fl;				/* Stream from fl_nm */
    char tok[TOK_LEN];			/* Token from parameter file,
					   indicates what comes next */

    /* Radar names in a combination */
    char nms[MAX_RADARS][PATH_LEN_MAX];

    int cvg;				/* Coverage flag, boolean */
    double ss;				/* Data weight multiplier */
    unsigned h;				/* Index in a coverage array */

    num_radars_max = num_radars_in;

    /* Initialize arrays to 1, i.e. all points are covered. */
    for (h = 0; h < CVG_HASH_SZ; h++) {
	cvg_opt_bg[h] = I_BAD;
	cvg_sub_bg[h] = I_BAD;
	cvg_opt_fil[h] = I_BAD;
	cvg_sub_fil[h] = I_BAD;
	cvg_bg[h] = I_BAD;
	cvg_fil[h] = I_BAD;
	sseq_trip[h] = NAN;
    }

    if ( !(fl = fopen(fl_nm, "r")) ) {
	fprintf(stderr, "Could not open calculation parameter file %s for "
		"reading.\n", fl_nm);
	return 0;
    }
    while (fscanf(fl, " " TOK_LEN_FMT, tok) == 1) {
	int r;				/* Radar index */
	int n;				/* Number of radars in combination */

	if ( strcmp(tok, "cvg_opt_bg:") == 0 ) {
	    if ( fscanf(fl, " %d", &n) != 1 ) {
		fprintf(stderr, "Failed to read radar count for cvg_opt_bg "
			"from %s\n", fl_nm);
		return 0;
	    }
	    for (r = 0 ; r < n; r++) {
		if ( fscanf(fl, " " PATH_LEN_FMT, nms[r]) != 1 ) {
		    fprintf(stderr, "Failed to read radar %d for cvg_opt_bg "
			    "from %s\n", r, fl_nm);
		    return 0;
		}
	    }
	    if ( fscanf(fl, " %d", &cvg) != 1 ) {
		fprintf(stderr, "Failed to read coverage value for for "
			"cvg_opt_bg from %s\n", fl_nm);
		return 0;
	    }
	    set_cvg_opt_bg(cvg, n, nms);
	} else if ( strcmp(tok, "cvg_sub_bg:") == 0 ) {
	    if ( fscanf(fl, " %d", &n) != 1 ) {
		fprintf(stderr, "Failed to read radar count for cvg_sub_bg "
			"from %s\n", fl_nm);
		return 0;
	    }
	    for (r = 0 ; r < n; r++) {
		if ( fscanf(fl, " " PATH_LEN_FMT, nms[r]) != 1 ) {
		    fprintf(stderr, "Failed to read radar %d for cvg_sub_bg "
			    "from %s\n", r, fl_nm);
		    return 0;
		}
	    }
	    if ( fscanf(fl, " %d", &cvg) != 1 ) {
		fprintf(stderr, "Failed to read coverage value for for "
			"cvg_sub_bg from %s\n", fl_nm);
		return 0;
	    }
	    set_cvg_sub_bg(cvg, n, nms);
	} else if ( strcmp(tok, "cvg_opt_fil:") == 0 ) {
	    if ( fscanf(fl, " %d", &n) != 1 ) {
		fprintf(stderr, "Failed to read radar count for cvg_opt_fil "
			"from %s\n", fl_nm);
		return 0;
	    }
	    for (r = 0 ; r < n; r++) {
		if ( fscanf(fl, " " PATH_LEN_FMT, nms[r]) != 1 ) {
		    fprintf(stderr, "Failed to read radar %d for cvg_opt_fil "
			    "from %s\n", r, fl_nm);
		    return 0;
		}
	    }
	    if ( fscanf(fl, " %d", &cvg) != 1 ) {
		fprintf(stderr, "Failed to read coverage value for for "
			"cvg_opt_fil from %s\n", fl_nm);
		return 0;
	    }
	    set_cvg_opt_fil(cvg, n, nms);
	} else if ( strcmp(tok, "cvg_sub_fil:") == 0 ) {
	    if ( fscanf(fl, " %d", &n) != 1 ) {
		fprintf(stderr, "Failed to read radar count for cvg_sub_fil "
			"from %s\n", fl_nm);
		return 0;
	    }
	    for (r = 0 ; r < n; r++) {
		if ( fscanf(fl, " " PATH_LEN_FMT, nms[r]) != 1 ) {
		    fprintf(stderr, "Failed to read radar %d for cvg_sub_fil "
			    "from %s\n", r, fl_nm);
		    return 0;
		}
	    }
	    if ( fscanf(fl, " %d", &cvg) != 1 ) {
		fprintf(stderr, "Failed to read coverage value for for "
			"cvg_sub_fil from %s\n", fl_nm);
		return 0;
	    }
	    set_cvg_sub_fil(cvg, n, nms);
	} else if ( strcmp(tok, "cvg_bg:") == 0 ) {
	    if ( fscanf(fl, " " PATH_LEN_FMT " %d" , nms[0], &cvg) != 2 ) {
		fprintf(stderr, "Could not read cvg_bg.\n");
		return 0;
	    }
	    set_cvg_bg(cvg, nms);
	} else if ( strcmp(tok, "cvg_fil:") == 0 ) {
	    if ( fscanf(fl, " " PATH_LEN_FMT " %d", nms[0], &cvg) != 2 ) {
		fprintf(stderr, "Could not read cvg_fil.\n");
		return 0;
	    }
	    set_cvg_fil(cvg, nms);
	} else if ( strcmp(tok, "sseq_trip:") == 0 ) {
	    if ( fscanf(fl, " %s %lf", nms[0], &ss) != 2 ) {
		fprintf(stderr, "Could not read sseq_trip.\n");
		return 0;
	    }
	    set_sseq_trip(ss, nms);
	} else if ( strcmp(tok, "#") == 0 ) {
	    int c;
	    for (c = fgetc(fl); (char)c != '\n' && c != EOF; c = fgetc(fl)) {
	    }
	}
    }

    return 1;
}

/*
   Set coverage arrays to default if not setting from coverage file.
 */

void Coverage_UseDefault(void)
{
    unsigned h;				/* Index in a coverage array */

    /* Initialize arrays to 1, i.e. all points are covered. */
    for (h = 0; h < CVG_HASH_SZ; h++) {
	cvg_opt_bg[h] = 1;
	cvg_sub_bg[h] = 1;
	cvg_opt_fil[h] = 1;
	cvg_sub_fil[h] = 1;
	cvg_bg[h] = 1;
	cvg_fil[h] = 1;
	sseq_trip[h] = 1.0;
    }
}

/*
   Compute hash for combination of num_radars radars with names nms.
   Combination must be available, or process exits.
 */

static unsigned int cvg_idx_s(int num_radars, char nms[MAX_RADARS][PATH_LEN_MAX])
{
    int n;				/* Radar loop index */
    int i;				/* Global radar index */
    unsigned h;				/* Index from a coverage array */

    for (n = 0, h = 0; n < num_radars; n++) {
	if ( (i = DDA_Get_Radar_Idx(nms[n])) == -1 ) {
	    fprintf(stderr, "Could not find storage location for %s. "
		    "Please check radar names in calc params file.\n", nms[n]);
	    exit(EXIT_FAILURE);
	}
	h |= 1 << i;
    }
    if ( h < 0 || h >= CVG_HASH_SZ ) {
	fprintf(stderr, "Could not find coverage value for %s. "
		"There is a bug in the program.\n", nms[n]);
	exit(EXIT_FAILURE);
    }
    return h;
}

/*
   Compute hash for combination of num_radars radars with indeces in idcs.
   Returns -1 on failure.
 */

static int cvg_idx_i(int num_radars, int *idcs)
{
    int n;				/* Radar index */
    unsigned h;				/* Index from a coverage array */

    for (n = 0, h = 0; n < num_radars; n++) {
	h |= 1 << idcs[n];
    }
    if ( h < 0 || h >= CVG_HASH_SZ ) {
	fprintf(stderr, "Could not find coverage value for radar combination");
	for (n = 0; n < num_radars; n++) {
	    fprintf(stderr, " %s", DDA_Get_Radar_Nm(idcs[n])
		    ? DDA_Get_Radar_Nm(idcs[n]) :  "unknown");
	}
	fprintf(stderr, "\n");
	return -1;
    }
    return h;
}

/*
   Set background coverage flag for a set of radars with optimal beam crossing
   angle. The group of radars denoted by the first num_radars names in nms will
   receive a coverage value of cvg. Returns 1/0 for success/failure.
 */

static void set_cvg_opt_bg(int cvg, int num_radars,
	char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(num_radars, nms);
    cvg_opt_bg[h] = cvg;
}

/*
   Set background coverage flag for a set of radars with sub-optimal beam
   crossing angle.
 */

static void set_cvg_sub_bg(int cvg, int num_radars,
	char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(num_radars, nms);
    cvg_sub_bg[h] = cvg;
}

/*
   Set filtered coverage flag for a set of radars with optimal beam
   crossing angle.
 */

static void set_cvg_opt_fil(int cvg, int num_radars,
	char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(num_radars, nms);
    cvg_opt_fil[h] = cvg;
}

/*
   Set filtered coverage flag for a set of radars with sub-optimal beam
   crossing angle.
 */

static void set_cvg_sub_fil(int cvg, int num_radars,
	char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(num_radars, nms);
    cvg_sub_fil[h] = cvg;
}

/* Set background coverage flag for one radar */ 
static void set_cvg_bg(int cvg, char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(1, nms);
    cvg_bg[h] = cvg;
}

/* Set filtered coverage flag for one radar */ 
static void set_cvg_fil(int cvg, char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(1, nms);
    cvg_fil[h] = cvg;
}

/* Set data weight multiplier */
static void set_sseq_trip(double ss, char nms[MAX_RADARS][PATH_LEN_MAX])
{
    unsigned h = cvg_idx_s(1, nms);
    sseq_trip[h] = ss;
}

/*
   Get background coverage flag for a set of radars with optimal beam crossing
   angle. The group of radars is denoted by the first num_radars indeces in
   idcs.
 */

int Get_Cvg_Opt_BG(int num_radars, int *idcs)
{
    unsigned h;

    if ( (h = cvg_idx_i(num_radars, idcs)) == -1 || cvg_opt_bg[h] == I_BAD ) {
	fprintf(stderr, "Opt BG coverage not set for");
	for (int n = 0; n < num_radars; n++) {
	    fprintf(stderr, " %s", DDA_Get_Radar_Nm(idcs[n])
		    ? DDA_Get_Radar_Nm(idcs[n]) :  "unknown");
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
    }
    return cvg_opt_bg[h];
}

/*
   Get background coverage flag for a set of radars with sub-optimal beam
   crossing angle. The group of radars is denoted by the first num_radars
   indeces in idcs.
 */

int Get_Cvg_Sub_BG(int num_radars, int *idcs)
{
    unsigned h;

    if ( (h = cvg_idx_i(num_radars, idcs)) == -1 || cvg_sub_bg[h] == I_BAD ) {
	fprintf(stderr, "Sub BG coverage not set for");
	for (int n = 0; n < num_radars; n++) {
	    fprintf(stderr, " %s", DDA_Get_Radar_Nm(idcs[n])
		    ? DDA_Get_Radar_Nm(idcs[n]) :  "unknown");
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
    }
    return cvg_sub_bg[h];
}

/*
   Get filtered coverage flag for a set of radars with optimal beam
   crossing angle.
 */

int Get_Cvg_Opt_Fil(int num_radars, int *idcs)
{
    unsigned h;

    if ( (h = cvg_idx_i(num_radars, idcs)) == -1 || cvg_opt_fil[h] == I_BAD ) {
	fprintf(stderr, "Opt fil coverage not set for");
	for (int n = 0; n < num_radars; n++) {
	    fprintf(stderr, " %s", DDA_Get_Radar_Nm(idcs[n])
		    ? DDA_Get_Radar_Nm(idcs[n]) :  "unknown");
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
    }
    return cvg_opt_fil[h];
}

/*
   Get filtered coverage flag for a set of radars with sub-optimal beam
   crossing angle.
 */

int Get_Cvg_Sub_Fil(int num_radars, int *idcs)
{
    unsigned h;

    if ( (h = cvg_idx_i(num_radars, idcs)) == -1 || cvg_sub_fil[h] == I_BAD ) {
	fprintf(stderr, "Sub fil coverage not set for");
	for (int n = 0; n < num_radars; n++) {
	    fprintf(stderr, " %s", DDA_Get_Radar_Nm(idcs[n])
		    ? DDA_Get_Radar_Nm(idcs[n]) :  "unknown");
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
    }
    return cvg_sub_fil[h];
}

/* Get background coverage flag radar r */ 
int Get_Cvg_BG(int r)
{
    unsigned h;

    if ( (h = cvg_idx_i(1, &r)) == -1 || cvg_bg[h] == I_BAD ) {
	fprintf(stderr, "BG coverage not set for %s.\n",
		DDA_Get_Radar_Nm(r) ? DDA_Get_Radar_Nm(r) :  "unknown");
	exit(EXIT_FAILURE);
    }
    return cvg_bg[h];
}

/* Get filtered coverage flag for radar r */ 
int Get_Cvg_Fil(int r)
{
    unsigned h;

    if ( (h = cvg_idx_i(1, &r)) == -1 || cvg_fil[h] == I_BAD ) {
	fprintf(stderr, "Fil coverage not set for %s.\n",
		DDA_Get_Radar_Nm(r) ? DDA_Get_Radar_Nm(r) :  "unknown");
	exit(EXIT_FAILURE);
    }
    return cvg_fil[h];
}

/* Get weight muliplier for radar r */ 
double Get_SSeq_Trip(int r)
{
    unsigned h;

    if ( (h = cvg_idx_i(1, &r)) == -1 || isnan(sseq_trip[h]) ) {
	fprintf(stderr, "SSeq coverage not set for %s.\n",
		DDA_Get_Radar_Nm(r) ? DDA_Get_Radar_Nm(r) :  "unknown");
	exit(EXIT_FAILURE);
    }
    return sseq_trip[h];
}

