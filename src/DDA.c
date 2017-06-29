/*
   Variational multiple-Doppler wind retrieval technique
   Developed by Alan Shapiro and Corey Potvin at Univ. of Oklahoma (OU);
   Coded by Corey Potvin (corey.potvin@noaa.gov)

   Code developed with funding from NSF grant ATM-0532107
   .	Title: Investigation of tornadic storms using Doppler polarimetric radar
   .	PIs:   A.V. Ryzhkov, T.-Y. Yu, and A. Shapiro

   Publications using this code must cite Potvin et al. (2012, JTECH) and,
   if the vorticity constraint was used, Shapiro et al. (2009, JTECH) as well
   (though may be better to cite both pubs even if vorticity constraint not
   used).

   For education or research use only.

   Search for "user-defined parameters" and "user-defined variables"
   Users may need to substantially modify other parts of code to suit their
   needs, especially subroutines for reading radar data files. Please send
   suggestions and bug reports to Corey.

   See README for more details.

   $Revision: 1.30 $ $Date: 2015/01/12 23:59:00 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <netcdf.h>
#include <unistd.h>
#include <setjmp.h>
#include "DDA.h"
#include "alloc.h"
#include "coverage.h"
#include "nnetcdf.h"
#include "geog_lib.h"
#include "geog_proj.h"

#define ERRCODE 2 // For netcdf error messages
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define PI 3.14159265358979323
#define Rad(X) (X*PI/180.0)
#define Deg(X) (X*180.0/PI)

/* International standard nautical mile = 1852 m = 1/60 deg of latitude */
#define NMILE 1852.0

// For Numerical Recipes conjugate gradient minimization code

#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define TOL 2.0e-4
#define NR_END 1
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FREE_ARG char*
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

int itmax_dbrent = 200;			/* Maximum number of iterations in
					   dbrent */
int itmax_frprmn0 = 200;		/* Number of iterations in frprmn
					   before the filtering round */
int itmax_frprmn1 = 10;			/* Number of iterations in frprmn
					   after the filtering round */


static int MM;				/* Total number of observation points,
					   = number of radars times number of
					   grid points actually used.*/

static double beta=1.0;
static double UT = 0, VT = 0;		/* Prescribed storm motion
					   (spatiotemporally constant) */
static double ss=1.0, ss2=1.0;
static int anel = 1;			/* For mass conservation constraint
					   (CalcDiv):
					   0 = Boussinesq approximation;
					   1 = anelastic approximation */
static int ncom;
static int laplace = 1;			/* For smoothness constraint (CalcSmooth):
					   0 = first-order derivatives;
					   1 = second-order derivatives */

static int weak_height=-1;		/* The index of the height equal to and
					   below that the sounding constraint is
					   weakened inside regions with greater
					   than 10 dbz. Set to -1 to prevent
					   implementation*/
static int upper_bc=1;			/* If 1 => impose upper BC w = 0
					   If -1 => ignore */
static int output_error=1;		/* If = 1, output verification stats
					   after each iteration (see bottom of
					   CalcCost) */
static int docount=1;
static int Vr_error=0;
static double C2, C3, C4, C5, C6, C7, C1a, C2a, C3a, C4a, C5a, C6a, C7a, C8a;

static double meanobserror_old=0.0; /* Old Mean obs error */
static double norm_merr_old = 0.0; /* Old Normalized Mass Cont error */
static double mass_norm_thresh = 0.3;   /* Normalized mass residual threshold */
static double meanobserror;		/* normalized obs error */
static double smoothtotal;		/* Fraction of Cost function due to smoothness */
static double obstotal;		    /* Fraction of Cost function due to observations */
static double masstotal;		/* Fraction of Cost function due to mass continuity */
static double vorttotal;		/* Fraction of Cost function due to voritcity */
static double Jtotal;           /* Total Cost Function minus the background */
static double maxobserror;      /* Max absolute obs error */

/* FILTER OPTIONS */
static unsigned filt_freq = 10;		/* Iteration frequency for the filter */
static enum {
    NO_FILTER,
    LEISE,
    LO_PASS
} filt_type = NO_FILTER;		/* Filter type */
static double filt_alpha = 0.15;	/* Alpha for low pass filter */
static int    Lei_steps  = 1;		/* Number of steps for Leise filter */
static int    vary_weights = 1;		/* If 1, then allow the weights to
					   vary */

/* Constraint weighting coefficients */
static double C1b;			/* Data constraint */
static double C2b;			/* Mass conservation */
static double C3b;			/* Vorticity */
static double C4b;			/* Horizontal smoothing */
static double C5b;			/* Vertical smoothing */
static double C6b;
static double C7b;
static double C8b;			/* Sounding constraint */

static double sss;
static double C1c, C2c, C3c, C4c, C5c, C6c, C7c, C8c, C1d, C2d, C3d, C4d, C5d;
static double C6d, C7d, C8d;
static double rmsVr;
static double *pcom, *xicom;

static int *map;
static int ***nn;
static int ***nn2;
static int ***coverage_bg;
static int ***coverage_fil;
static int **obs_index;			/* Index that keeps track of the obs
					   index for each height for each
					   radar */
static float ***vort;
static float ***DzetaDt;
static float ***div2;
static float ***UUU;
static float ***VVV;
static float ***dudx;
static float ***dudy;
static float ***dvdx;
static float ***dvdy;
static float ***dwdz;
static float ***dwdz2;
static float ***dwdz3;
static float ***dwdx2;
static float ***dwdy2;
static float ***dudx2;
static float ***dudy2;
static float ***dvdx2;
static float ***dvdy2;
static float ***dvdz2;
static float ***dudz2;
static float ***dwdy;
static float ***dvdz;
static float ***dudz;
static float ***dwdx;
static double *rho;
static double ***maxdbz;
static double *W1a;
static float *SS;
static float *V_obs;
static float *refl_obs;
static float *term_vel;			/* Hydrometeor terminal veloctiy */
static double *p;
static double *xp;
static double *truth;
static double *p_old;
static double *ub;
static double *vb;

static int *obs_mask;
static float *xpos;
static float *ypos;
static float *zpos;
static double *Rmaxx;
static double *Rmaxx2;
static float *azims;
static double *ran;
static float *elevs;
static double *vr_retr;

static double CalcTermVel(double, double);
static void init_wind3(int, size_t, size_t, size_t, float *, float *, float *,
	float *, double *, double *, float *, float *, float, float,
	double);
static void ReadPyART (int, char [MAX_RADARS][PATH_LEN_MAX], const char *,
	const char *, double);
static void Filt2d(double **, double **, int **, int, int, float);
static void Velfill(double *, int, int, int, int, int, int, float);
static void VelLefilt(double *, int , int, int , int);
static void writeout(double *);
static double CalcCost(double[]);
static double CalcDiv(double *);
static double CalcVortTend(double *);
static double CalcSmooth(double *);
static double CalcBG(double *);
static void ReadBG(char *);
static void CalcGrad(double *, double *);
static double dbrent(double, double, double,
	double (*)(double),
	double (*f)(double), double, double *);
static void mnbrak(double *, double *, double *, double *, double *, double *,
	double (*)(double));
static void dlinmin(double *, double *, int, double (*)(double *),
	void (*)(double *, double *), double *);
static void frprmn(double *, int, double, int *, int, double (*)(double *),
	void (*)(double *, double *));
static double *vector(long int, long int);
static void free_vector(double *, long, long);
static void nrerror(const char *);
static double f1dim(double);
static double df1dim(double);
static void GradCheck(int, double[]);
static double (*nrfunc)(double *);
static void (*nrdfun)(double *, double *);

static double dist2(const double, const double, const double, const double);

/* Leise filter function, define in leise_filt.f */
void t5fltr_(float *, int *, int *, int *, int *);

static int *ialloc(size_t, const char *);
static float *falloc(size_t, const char *);
static double *dalloc(size_t, const char *);
static int **ialloc2(long, long, const char *);
static void ifree2(int **);
static double **dalloc2(long, long, const char *);
static void dfree2(double **);
static int ***ialloc3(long, long, long, const char *);
static float ***falloc3(long, long, long, const char *);
static double ***dalloc3(long, long, long, const char *);
static float ****calloc4f(long , long , long , long , float *, char *);
static void free4f(float ****);

/* Names of output files */
char frprmn_fl_nm[PATH_LEN_MAX];	/* Created by frprmn */
char writeout_fl_nm[PATH_LEN_MAX];	/* Created by writeout */

/* Grid parameters. Distances in meters, angles in degrees. */
static int num_x, num_y, num_z;		/* Number of x, y, z coordinates in
					   analysis domain */
static int num_dx, num_dy, num_dz;	/* Number of steps, one less than number
					   of coordinates */
static int N;				/* Total number of coordinates */
static double x_min, y_min, z_min;	/* Coordinates of west, south, top */
static double x_max, y_max, z_max;	/* Coordinates of east, north, bottom */
static double x_inc, y_inc, z_inc;	/* Coordinate spacing */
static float grid_lon, grid_lat;	/* Longitude,latitude at x = y = 0 */
static float grid_alt;			/* Altitude at x = y = 0 */
static double cutoff = 0.0;		/* Optional height below which to omit
					   observations from analysis
					   (e.g., data-denial experiments) */
/* Radar locations */
static int num_radar;			/* Number of radars */
static double radx[MAX_RADARS], rady[MAX_RADARS];
static float radz[MAX_RADARS];

/* Radar names = base names of files that provide radar data */
static char radar_nms[MAX_RADARS][PATH_LEN_MAX];

int main(int argc, char *argv[])
{
    char *argv0 = argv[0];		/* This command */
    char *param_fl_nm;			/* Parameter file name */
    FILE *param_fl;			/* Parameter file */
    char tok[TOK_LEN];			/* Token from parameter file,
					   indicates what comes next */
    char dir[PATH_LEN_MAX];		/* Execution directory */

    /* OPAWS files to read */
    char opaws_fl_nms[MAX_RADARS][PATH_LEN_MAX];

    /* Name of file with calculation parameters */
    char calc_param_fl_nm[PATH_LEN_MAX];

    char bg_fl[PATH_LEN_MAX];		/* Background file */
    FILE *calc_param_fl;		/* File with calculation parameters */
    double H = 10000.0;			/* Scale height for density */
    double rho_0 = 1.0;			/* Density at z = 0.0 */
    double min_CBA = 20.0;		/* dual-Doppler domain criteria:
					   minimum threshold on radar cross-beam
					   angle */
    double max_dist = 10.0;		/* dual-Doppler domain criteria:
					   there must be at least one obs from
					   each radar within this distance of
					   the analysis point */

    int m, n, k;
    char refl_nm[TOK_LEN];		/* Name of reflectivity variable in
					   input file. */
    char vt_nm[TOK_LEN];		/* Name of velocity variable in
					   input file. */
    FILE *W;
    int iter;
    int read_dataweights = 2;		/* For weights in data constraint:
					   0 = calculate and output to file
					   1 = read from file
					   2 = weight all obs equally */

    /* Initialize */
    printf("DDA %s\n", DDA_VERSION);
    if ( argc != 2 ) {
	fprintf(stderr, "Usage: %s parameter_file\n", argv0);
	exit(EXIT_FAILURE);
    }
    param_fl_nm = argv[1];
    x_min = x_inc = y_min = y_inc = z_min = z_inc = NAN;
    grid_lon = grid_lat = grid_alt = NAN;
    strcpy(dir, ".");
    strcpy(refl_nm, "DZ");
    strcpy(vt_nm, "VT");
    strcpy(bg_fl, "");
    strcpy(calc_param_fl_nm, "");
    strcpy(frprmn_fl_nm, "");
    strcpy(writeout_fl_nm, "");

    /*
       Initialize constraint weighting coefficients with bogus values.
       These will be set from calc_params file.
       If not, they will be set to proper defaults below.
     */

    C1b = C2b = C3b = C4b = C7b = C5b = C8b = NAN;

    /*
       Read parameter file.
       The parameter file must contain the following sets of tokens and
       associated values:

       |	x: x_min x_inc num_x
       |	y: y_min y_inc num_y
       |	z: z_min z_inc num_z
       |	grid:	lon lat
       |	center:	lon lat
       |	dir: directory_path
       |	opaws: num_files file_name file_name file_name ...
       |	refl: reflectivity_variable
       |	bgfile: background_file
       |	writeout: writeout_file
       |	frprmn_out: frprmn_file
       |	min_cba: min_beam_cross
       |	# comment ...

       Lines can be in any order. Tokens and values can be separated by
       any combination of white space.

       If "#" is encountered, input is skipped to next newline.

       dir specifies working directory for process execution. Default is
       current working directory.

       There must be one opaws file for each radar.
     */

    if ( !(param_fl = fopen(param_fl_nm, "r")) ) {
	fprintf(stderr, "%s: could not open parameter file %s for reading.\n",
		argv0, param_fl_nm);
	exit(EXIT_FAILURE);
    }
    while (fscanf(param_fl, " " TOK_LEN_FMT, tok) == 1) {
	if ( strcmp(tok, "x:") == 0 ) {
	    if ( fscanf(param_fl, " %lf %lf %d",
			&x_min, &x_inc, &num_x) != 3 ) {
		fprintf(stderr, "%s: failed to read x grid specificiation "
			"from %s\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "y:") == 0 ) {
	    if ( fscanf(param_fl, " %lf %lf %d",
			&y_min, &y_inc, &num_y) != 3 ) {
		fprintf(stderr, "%s: failed to read y grid specificiation "
			"from %s\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "z:") == 0 ) {
	    if ( fscanf(param_fl, " %lf %lf %d",
			&z_min, &z_inc, &num_z) != 3 ) {
		fprintf(stderr, "%s: failed to read z grid specificiation "
			"from %s\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "grid:") == 0 ) {
	    if ( fscanf(param_fl, " %f %f %f",
			&grid_lon, &grid_lat, &grid_alt) != 3 ) {
		fprintf(stderr, "%s: failed to read grid origin from %s\n",
			argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	    grid_lon = Rad(grid_lon);
	    grid_lat = Rad(grid_lat);
	} else if ( strcmp(tok, "dir:") == 0 ) {
	    if ( fscanf(param_fl, PATH_LEN_FMT, dir) != 1 ) {
		fprintf(stderr, "%s: could not read execution directory "
			"from %s.\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "opaws:") == 0 ) {
	    if ( fscanf(param_fl, " %d", &num_radar) != 1 ) {
		fprintf(stderr, "%s: failed to read z grid specificiation "
			"from %s\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	    for (n = 0; n < num_radar; n++) {
		if ( fscanf(param_fl, PATH_LEN_FMT " " PATH_LEN_FMT,
			    opaws_fl_nms[n], radar_nms[n]) != 2 ) {
		    fprintf(stderr, "%s: failed to read name of OPAWS file and "
			    "radar name index %d from parameter file %s.\n",
			    argv0, n, param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    }
	} else if ( strcmp(tok, "refl:") == 0 ) {
	    if ( fscanf(param_fl, TOK_LEN_FMT, refl_nm) != 1 ) {
		fprintf(stderr, "%s: failed to read name of reflectivity "
			"variable from parameter file %s.\n",
			argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "vt:") == 0 ) {
	    if ( fscanf(param_fl, TOK_LEN_FMT, vt_nm) != 1 ) {
		fprintf(stderr, "%s: failed to read name of velocity "
			"variable from parameter file %s.\n",
			argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "bgfile:") == 0 ) {
	    if ( fscanf(param_fl, PATH_LEN_FMT, bg_fl) != 1 ) {
		fprintf(stderr, "%s: failed to read name of background file"
			"from parameter file %s.\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "writeout:") == 0 ) {
	    if ( fscanf(param_fl, PATH_LEN_FMT, writeout_fl_nm) != 1 ) {
		fprintf(stderr, "%s: failed to read name of writeout output "
			"file from parameter file %s.\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "frprmn_out:") == 0 ) {
	    if ( fscanf(param_fl, PATH_LEN_FMT, frprmn_fl_nm) != 1 ) {
		fprintf(stderr, "%s: failed to read name of frprmn output "
			"file from parameter file %s.\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "min_cba:") == 0 ) {
	    if ( fscanf(param_fl, " %lf", &min_CBA) != 1 ) {
		fprintf(stderr, "%s: failed to read minimum beam crossing"
			" angle from %s\n", argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	    if ( min_CBA < 0.0 || min_CBA > 90.0 ) {
		fprintf(stderr, "%s: Minimum beam crossing angle must be"
			" between 0 and 90 degrees, not %f.\n", argv0, min_CBA);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "calc_params:") == 0 ) {
	    if ( fscanf(param_fl, PATH_LEN_FMT, calc_param_fl_nm) != 1 ) {
		fprintf(stderr, "%s: failed to read name of calculation "
			"parameter file from parameter file %s.\n",
			argv0, param_fl_nm);
		exit(EXIT_FAILURE);
	    }
	} else if ( strcmp(tok, "#") == 0 ) {
	    for (int c = fgetc(param_fl);
		    (char)c != '\n' && c != EOF;
		    c = fgetc(param_fl)) {
	    }
	}
    }
    fclose(param_fl);

    /* Validate input from parameter file */
    if ( !isfinite(x_min + x_inc) || x_inc <= 0.0 || num_x == 0 ) {
	fprintf(stderr, "%s: could not find x grid specification in %s.\n",
		argv0, param_fl_nm);
	exit(EXIT_FAILURE);
    }
    if ( !isfinite(y_min + y_inc) || y_inc <= 0.0 || num_y == 0 ) {
	fprintf(stderr, "%s: could not find y grid specification in %s.\n",
		argv0, param_fl_nm);
	exit(EXIT_FAILURE);
    }
    if ( !isfinite(z_min + z_inc) || z_inc <= 0.0 || num_z == 0 ) {
	fprintf(stderr, "%s: could not find z grid specification in %s.\n",
		argv0, param_fl_nm);
	exit(EXIT_FAILURE);
    }
    if ( !isfinite(grid_lon + grid_lat + grid_alt) ) {
	fprintf(stderr, "%s: could not find grid origin in %s.\n",
		argv0, param_fl_nm);
	exit(EXIT_FAILURE);
    }
    if ( fabs(grid_lat) > 90.0 ) {
	fprintf(stderr, "%s: grid latitude %f > 90.0 degrees. "
		"Are longitude and latitude transposed.", argv0, grid_lat);
	exit(EXIT_FAILURE);
    }
    if ( strlen(frprmn_fl_nm) == 0 ) {
	fprintf(stderr, "%s: no name given for output file in frprmn.\n",
		argv0);
	exit(EXIT_FAILURE);
    }
    if ( strlen(writeout_fl_nm) == 0 ) {
	fprintf(stderr, "%s: no name given for output file in writeout "
		"function.\n", argv0);
	exit(EXIT_FAILURE);
    }
    if ( num_radar <= 0 ) {
	fprintf(stderr, "%s: number of radars must be positive, "
		"not %d\n", argv0, num_radar);
	exit(EXIT_FAILURE);
    }
    for (n = 0; n < num_radar; n++) {
	if ( !opaws_fl_nms[n] ) {
	    fprintf(stderr, "%s: no name for radar index %d\n",
		    argv0, n);
	    exit(EXIT_FAILURE);
	}
    }

    /*
       If specified, read calculation parameters from calc_params file.
       See sample file calc_params.dda for explanation.
       The parameter file may contain the following sets of tokens and
       associated values:

       |	anel: flag
       |	laplace: flag
       |	read_dataweights: flag
       |	max_dist: float
       |	cutoff: float
       |	UT: float
       |	VT: float
       |	output_error: flag
       |	weak_height: integer
       |	upper_bc: flag
       |	itmax_frprmn: integer integer
       |	itmax_dbrent: integer
       |	C1b: float
       |	C2b: float
       |	C3b: float
       |	C4b: float
       |	C5b: float
       |	C8b: float
       |	# comment ...

       If a token is absent, default will prevail.
       Lines can be in any order. Tokens and values can be separated by
       any combination of white space.

       If "#" is encountered, input is skipped to next newline.
     */

    if ( strlen(calc_param_fl_nm) > 0 ) {
	if ( !(calc_param_fl = fopen(calc_param_fl_nm, "r")) ) {
	    fprintf(stderr, "%s: could not open calc params file %s for "
		    "reading.\n", argv0, param_fl_nm);
	    exit(EXIT_FAILURE);
	}
	printf("%s: reading calculation parameters from %s.\n",
		argv0, calc_param_fl_nm);
	while (fscanf(calc_param_fl, " %s", tok) == 1) {
	    if ( strcmp(tok, "anel:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &anel) != 1 ) {
		    fprintf(stderr, "%s: failed to read flag for mass "
			    "conservation constraint from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "laplace:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &laplace) != 1 ) {
		    fprintf(stderr, "%s: failed to read flag for smoothness "
			    "constraint from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "read_dataweights:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &read_dataweights) != 1 ) {
		    fprintf(stderr, "%s: failed to read flag indicating how "
			    "to get data weights from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "max_dist:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &max_dist) != 1 ) {
		    fprintf(stderr, "%s: failed to read maximum distance "
			    "between radar coverage point and grid point "
			    "from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "cutoff:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &cutoff) != 1 ) {
		    fprintf(stderr, "%s: failed to read cutoff height "
			    "from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "UT:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &UT) != 1 ) {
		    fprintf(stderr, "%s: failed to read u component of "
			    "storm motion from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "VT:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &VT) != 1 ) {
		    fprintf(stderr, "%s: failed to read v component of "
			    "storm motion from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "output_error:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &output_error) != 1 ) {
		    fprintf(stderr, "%s: failed to read flag requesting "
			    "verification statistics after each iteration "
			    "from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "weak_height:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &weak_height) != 1 ) {
		    fprintf(stderr, "%s: failed to read index of weak "
			    "height level from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "upper_bc:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &upper_bc) != 1 ) {
		    fprintf(stderr, "%s: failed to read flag requesting "
			    "upper boundary condition from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "itmax_frprmn:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d %d", &itmax_frprmn0,
			    &itmax_frprmn1 ) != 2 ) {
		    fprintf(stderr, "%s: failed to read maximum iteration "
			    "counts for frprmn from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "itmax_dbrent:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &itmax_dbrent) != 1 ) {
		    fprintf(stderr, "%s: failed to read maximum iteration "
			    "count for dbrent from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C1b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C1b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C1b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C2b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C2b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C2b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C3b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C3b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C3b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C4b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C4b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C4b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C5b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C5b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C5b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "C8b:") == 0 ) {
		if ( fscanf(calc_param_fl, " %lf", &C8b) != 1 ) {
		    fprintf(stderr, "%s: failed to read C8b from %s\n",
			    argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "filter:") == 0 ) {
		char filt_type_s[TOK_LEN];
		int nin;		/* Return value from fscanf */

		nin = fscanf(calc_param_fl, " %d", &filt_freq);
		if ( nin == 1 ) {
		    /*
		       Have filter frequency.
		       Get filter type and parameters.
		     */

		    nin = fscanf(calc_param_fl, " " TOK_LEN_FMT, filt_type_s);
		    if ( nin == 1 ) {
			if ( strcmp(filt_type_s, "Leise") == 0 ) {
			    filt_type = LEISE;
			    nin = fscanf(calc_param_fl, " %d", &Lei_steps);
			    if ( nin != 1 ) {
				fprintf(stderr, "%s: could not read Leise "
					"filter step count from %s\n",
					argv0, calc_param_fl_nm);
				exit(EXIT_FAILURE);
			    }
			} else if ( strcmp(filt_type_s, "low-pass") == 0 ) {
			    filt_type = LO_PASS;
			    nin = fscanf(calc_param_fl, " %lf", &filt_alpha);
			    if ( nin != 1 )  {
				fprintf(stderr, "%s: could not read low pass "
					"filter alpha from %s\n",
					argv0, calc_param_fl_nm);
				exit(EXIT_FAILURE);
			    }
			} else {
			    fprintf(stderr, "%s: Unknown filter type %s. "
				    "Filter type must be \"none\", \"Leise\", "
				    "or \"low-pass\"\n", argv0, filt_type_s);
			    exit(EXIT_FAILURE);
			}
		    }
		} else {
		    /*
		       No filter frequency given. Check for no filter.
		     */

		    nin = fscanf(calc_param_fl, " " TOK_LEN_FMT,
			    filt_type_s);
		    if ( nin == 1 && strcmp(filt_type_s, "none") == 0 ) {
			filt_type = NO_FILTER;
		    } else {
			fprintf(stderr, "%s: failed to read filter "
				    "parameters from %s\n",
				    argv0, calc_param_fl_nm);
			    exit(EXIT_FAILURE);
		    }
		}
	    } else if ( strcmp(tok, "vary_weights:") == 0 ) {
		if ( fscanf(calc_param_fl, " %d", &vary_weights) != 1 ) {
		    fprintf(stderr, "%s: failed to read vary_weights "
			    "parameter from %s\n", argv0, calc_param_fl_nm);
		    exit(EXIT_FAILURE);
		}
	    } else if ( strcmp(tok, "#") == 0 ) {
		for (int c = fgetc(param_fl);
			(char)c != '\n' && c != EOF;
			c = fgetc(param_fl)) {
		}
	    }
	}
	fclose(calc_param_fl);
	if ( !Coverage_Read(calc_param_fl_nm, num_radar) ) {
	    fprintf(stderr, "%s: failed to read coverage parameters from %s\n",
		    argv0, calc_param_fl_nm);
	    exit(EXIT_FAILURE);
	}
    } else {
	Coverage_UseDefault();
    }
    if ( anel != 0 && anel != 1 ) {
	fprintf(stderr, "%s: anel flag must be 0 or 1\n", argv0);
	exit(EXIT_FAILURE);
    }
    if ( laplace != 0 && laplace != 1 ) {
	fprintf(stderr, "%s: laplace flag must be 0 or 1\n", argv0);
	exit(EXIT_FAILURE);
    }
    if ( read_dataweights != 0 && read_dataweights != 1
	    && read_dataweights != 2 ) {
	fprintf(stderr, "%s: read dataweights flag must be one of 0, 1, or 2\n",
		argv0);
	exit(EXIT_FAILURE);
    }
    if ( output_error != 0 && output_error != 1 ) {
	fprintf(stderr, "%s: output error flag must be 0 or 1\n", argv0);
	exit(EXIT_FAILURE);
    }
    if ( upper_bc != -1 && upper_bc != 1 ) {
	fprintf(stderr, "%s: upper boundary condition flag must be 0 or 1\n",
		argv0);
	exit(EXIT_FAILURE);
    }
    int have_filter = 0;
    switch (filt_type) {
	case NO_FILTER:
	    have_filter = 1;
	    break;
	case LEISE:
	    have_filter = 1;
	    break;
	case LO_PASS:
	    have_filter = 1;
	    break;
    }
    if ( !have_filter ) {
	fprintf(stderr, "%s: fitler type is undefined.\n", argv0);
	exit(EXIT_FAILURE);
    }
    if ( Lei_steps < 0 ) {
	fprintf(stderr, "%s: number of steps in Leise filter cannot be "
		"negative", argv0);
	exit(EXIT_FAILURE);
    }

    if ( strcmp(dir, ".") != 0 ) {
	printf("Changing working directory to %s\n", dir);
	if ( chdir(dir) == -1 ) {
	    fprintf(stderr, "%s: could not change working directory to %s\n",
		    argv0, dir);
	    exit(EXIT_FAILURE);
	}
    }

    /* Analysis domain limits (domain origin specified further down) */
    num_dx = num_x - 1;
    num_dy = num_y - 1;
    num_dz = num_z - 1;
    x_max = x_min + num_dx * x_inc;
    y_max = y_min + num_dy * y_inc;
    z_max = z_min + num_dz * z_inc;
    N = num_x * num_y * num_z;

    /* Allocate arrays */
    nn = ialloc3(num_x, num_y, num_z, "nn");
    nn2 = ialloc3(num_x, num_y, num_z, "nn2");
    obs_index = ialloc2(num_z, MAX_RADARS, "obs_index");
    vort = falloc3(num_x, num_y, num_z, "vort");
    DzetaDt = falloc3(num_x, num_y, num_z, "DzetaDt");
    div2 = falloc3(num_x, num_y, num_z, "div2");
    UUU = falloc3(num_x, num_y, num_z, "UUU");
    VVV = falloc3(num_x, num_y, num_z, "VVV");
    dudx = falloc3(num_x, num_y, num_z, "dudx");
    dudy = falloc3(num_x, num_y, num_z, "dudy");
    dvdx = falloc3(num_x, num_y, num_z, "dvdx");
    dvdy = falloc3(num_x, num_y, num_z, "dvdy");
    dwdz = falloc3(num_x, num_y, num_z, "dwdz");
    dwdz2 = falloc3(num_x, num_y, num_z, "dwdz2");
    dwdz3 = falloc3(num_x, num_y, num_z, "dwdz3");
    dwdx2 = falloc3(num_x, num_y, num_z, "dwdx2");
    dwdy2 = falloc3(num_x, num_y, num_z, "dwdy2");
    dudx2 = falloc3(num_x, num_y, num_z, "dudx2");
    dudy2 = falloc3(num_x, num_y, num_z, "dudy2");
    dvdx2 = falloc3(num_x, num_y, num_z, "dvdx2");
    dvdy2 = falloc3(num_x, num_y, num_z, "dvdy2");
    dvdz2 = falloc3(num_x, num_y, num_z, "dvdz2");
    dudz2 = falloc3(num_x, num_y, num_z, "dudz2");
    dwdy = falloc3(num_x, num_y, num_z, "dwdy");
    dvdz = falloc3(num_x, num_y, num_z, "dvdz");
    dudz = falloc3(num_x, num_y, num_z, "dudz");
    dwdx = falloc3(num_x, num_y, num_z, "dwdx");
    rho = dalloc(1000, "rho");
    maxdbz = dalloc3(num_x, num_y, num_z, "maxdbz");
    xp = dalloc(3 * N + num_z * num_y + 2, "xp");
    p_old = dalloc(3 * N + 1, "p_old");
    ub = dalloc(N + 1, "ub");
    vb = dalloc(N + 1, "vb");
    truth = dalloc(3 * N + 1, "truth");
    coverage_bg = ialloc3(num_dx + 3, num_dy + 3, num_dz + 1, "coverage_bg");
    coverage_fil = ialloc3(num_dx + 3, num_dy + 3, num_dz + 1, "coverage_fil");

    printf("x: %lf to %lf in %d steps of %lf\n", x_min, x_max, num_dx, x_inc);
    printf("y: %lf to %lf in %d steps of %lf\n", y_min, y_max, num_dy, y_inc);
    printf("z: %lf to %lf in %d steps of %lf\n", z_min, z_max, num_dz, z_inc);
    printf("z_min=%g, cutoff=%g \n", z_min, cutoff);
    printf("UT=%g, VT=%g \n", UT, VT);

    if (anel==1) {
	printf("Anelastic mass cons\n");
    } else {
	printf("Boussinesq mass cons\n");
    }
    if (laplace==1) {
	printf("Second-order smoothness constraint\n");
    } else {
	printf("First-order smoothness constraint\n");
    }
    printf("Minimum beam crossing angle = %g degrees\n", min_CBA);

    /*
       Create array for mapping between 3-D and 1-D arrays (required since
       minimization routine frprmn takes in 1-D array argument).
     */

    n=0;
    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {
		n++;
		nn[i][j][k]=n;
	    }
	}
    }

    /* Read in radar data */
    ReadPyART(num_radar, opaws_fl_nms, refl_nm, vt_nm, min_CBA);
    /* ReadOPAWS(num_radar, opaws_fl_nms, refl_nm, vt_nm, min_CBA); */

    /* Allocate arrays whose size depends on number of observations read, MM */
    p = dalloc(3 * MM, "p");
    W1a = dalloc(MM, "W1a");
    obs_mask = ialloc(MM, "obs_mask");
    vr_retr = dalloc(MM, "vr_retr");
    Rmaxx = dalloc(MM, "Rmaxx");
    Rmaxx2 = dalloc(MM, "Rmaxx2");
    ran = dalloc(MM, "ran");

    /*
       Compute vertical density profile and derivative for anelastic
       continuity constraint.
     */

    for (k=0; k < num_z; k++) {
	double z = z_min + k * z_inc;
	rho[k] = rho_0 * exp(-z / H);
    }

    /* Set constraint weighting coefficients if not set from calc_params. */
    if ( !isfinite(C1b) ) {
	C1b = 1.0;
    }
    if ( !isfinite(C2b) ) {
	C2b = 5.0e4;
    }
    if ( !isfinite(C3b) ) {
	C3b = 0;
    }
    if ( !isfinite(C4b) ) {
	C4b = 5e6;
    }
    if ( !isfinite(C7b) ) {
	C7b = C4b;
    }
    if ( !isfinite(C5b) ) {
	C5b = 5e6*pow(z_inc / PI, 4) / N;
    }
    C6b = C5b;
    if ( !isfinite(C8b) ) {
	if ( strlen(bg_fl) > 0 ) {
	    C8b = 0.01;
	    ReadBG(bg_fl);
	} else {
	    C8b = 0.0;
	}
    } else if ( C8b == 0.0 ) {
	if ( strlen(bg_fl) > 0 ) {
	    fprintf(stderr, "%s: C8b cannot be 0.0 if using background file.\n",
		    argv0);
	}
    } else {
	if ( strlen(bg_fl) == 0 ) {
	    fprintf(stderr, "%s: C8b is not zero but no background file "
		    "specified.\n", argv0);
	}
	ReadBG(bg_fl);
    }

    /*
       Cost function computation; executed here to compute rmsVr
       [used to normalize data constraint weighting coefficient(s)].
     */

    CalcCost(xp);

    C1a=C1b; C2a=C2b; C3a=C3b; C4a=C4b; C5a=C5b; C6a=C6b; C7a=C7b; C8a=C8b;

    printf("C1b=%g C2b=%g C3b=%g C4b=%g C5b=%g C6b=%g C7b=%g C8b=%g\n",
	    C1b, C2b, C3b, C4b, C5b, C6b, C7b, C8b);
    printf("C2a=%g C3a=%g C4a=%g C5a=%g C6a=%g C7a=%g C8a=%g\n",
	    C2a, C3a, C4a, C5a, C6a, C7a, C8a);

    /* Weighting coefficients within buffer domain */
    C1c=0; C2c=C2b; C3c=C3b; C4c=C4b; C5c=C5b; C6c=C6b; C7c=C7b; C8c=C8b;
    C1d=0; C2d=C2c; C3d=C3c; C4d=C4c; C5d=C5c; C6d=C6c; C7d=C7c; C8d=C8c;

    /* Compute or read observational weights for data constraint */
    switch (read_dataweights) {
	case 0:
	    printf("Calculating data weights....");
	    if ( !(W = fopen("W1.dat", "w")) ) {
		printf("Could not open 'W1.dat' files.  Exiting program.\n");
		exit(EXIT_FAILURE);
	    }
	    for (m = 0; m < MM; m++) {
		fprintf(W, "%g\n", SS[m]);
		W1a[m] = C1b * SS[m] / (rmsVr * rmsVr);
	    }
	    fclose(W);
	    break;
	case 1:
	    printf("Reading in data weights....\n");
	    if ( !(W = fopen("W1.dat", "r")) ) {
		printf("Could not open 'W1.dat' files!  Exiting program.\n");
		exit(EXIT_FAILURE);
	    }
	    for (m=0; m<MM; m++) {
		if ( fscanf(W, " %g", SS + m) != 1 ) {
		    fprintf(stderr, "Could not read element %d of W1.dat", m);
		    exit(EXIT_FAILURE);
		}
		W1a[m] = C1b * SS[m] / (rmsVr * rmsVr);
	    }
	    fclose(W);
	    break;
	case 2:
	    printf("Weighting all obs equally! \n");
	    for (m=0; m<MM; m++) {
		W1a[m] = C1b * SS[m] / (rmsVr * rmsVr);
	    }
	    break;
    }

    /*
       Verify adjoint (CalcGrad) consistent with cost function (CalcCost)
       - want to see fr decrease to 1.0000 then zero
     */
    GradCheck(3*N,xp);

    /* Minimize J and output analysis */
    printf("C2b = %g C5b = %g %s\n", C2b, C5b,
	    vary_weights ? "Varying weights" : "Not varying weights.");
    printf("xp: %g %g %g %g %g %g\n",
	    xp[nn[36][34][7]], xp[nn[37][34][7]], xp[nn[38][34][7]],
	    xp[nn[39][34][7]], xp[nn[40][34][7]], xp[nn[41][34][7]]);
    printf("Starting minimization....\n");
    frprmn(xp, 3*N, TOL, &iter, itmax_frprmn0, CalcCost, CalcGrad);
    printf("Finished\n");
    printf("C2b = %g C5b = %g %s\n", C2b, C5b,
	    vary_weights ? "Varying weights" : "Not varying weights.");
    printf("xp: %g %g %g %g %g %g\n",
	    xp[nn[36][34][7]], xp[nn[37][34][7]], xp[nn[38][34][7]],
	    xp[nn[39][34][7]], xp[nn[40][34][7]], xp[nn[41][34][7]]);
    writeout(xp);
    printf("Filter applied\n");
    printf("C2b = %g C5b = %g %s\n", C2b, C5b,
	    vary_weights ? "Varying weights" : "Not varying weights.");
    printf("xp: %g %g %g %g %g %g\n",
	    xp[nn[36][34][7]], xp[nn[37][34][7]], xp[nn[38][34][7]],
	    xp[nn[39][34][7]], xp[nn[40][34][7]], xp[nn[41][34][7]]);

    /* Reset weight constraint weighting coefficients */
    C1b = 1;
    C2b = 1;
    C3b = 0;
    C4b = C7b = 5e3;
    C5b = C6b = 5e3 * pow(z_inc / PI, 4) / N;
    C8b = 0.0;

    /*
       Cost function computation; executed here to compute rmsVr
       [used to normalize data constraint weighting coefficient(s)].
     */

    CalcCost(xp);

    C1a = C1b;
    C2a = C2b;
    C3a = C3b;
    C4a = C4b;
    C5a = C5b;
    C6a = C6b;
    C7a = C7b;
    C8a = C8b;

    /* weighting coefficients within buffer domain */
    C1c = 0;
    C2c = C2b;
    C3c = C3b;
    C4c = C4b;
    C5c = C5b;
    C6c = C6b;
    C7c = C7b;
    C8c = C8b;
    C1d = 0;
    C2d = C2c;
    C3d = C3c;
    C4d = C4c;
    C5d = C5c;
    C6d = C6c;
    C7d = C7c;
    C8d = C8c;

    /*
       Verify adjoint (CalcGrad) consistent with cost function (CalcCost)
       - want to see fr decrease to 1.0000 then zero
     */

    GradCheck(3 * N, xp);
    frprmn(xp, 3 * N, TOL, &iter, itmax_frprmn1, CalcCost, CalcGrad);
    printf("Finished\n");
    printf("%g %g %g %g %g %g\n",
	    xp[nn[36][34][7]], xp[nn[37][34][7]], xp[nn[38][34][7]],
	    xp[nn[39][34][7]], xp[nn[40][34][7]], xp[nn[41][34][7]]);
}

/*
   Return name of radar r, or NULL if r is out of bounds. Caller should not
   modify return value
 */ 

char *DDA_Get_Radar_Nm(int r)
{
    if ( r < 0 || r >= num_radar ) {
	return NULL;
    }
    return radar_nms[r];
}

/*
   Return index of nm in global radar_nms array, or -1 is nm is not in
   radar_nms
 */

int DDA_Get_Radar_Idx(const char *nm)
{
    int r;				/* Index in radar_nms */

    for (r = 0; r < num_radar; r++) {
	if ( strcmp(nm, radar_nms[r]) == 0 ) {
	    return r;
	}
    }
    return -1;
}

void ReadBG(char *bg_fl)
{
    int i, j, k, n;
    double USND[num_z], VSND[num_z], z;
    int num_snd;			/* Number of sounding levels */
    int num_snd_max = 1000;
    float height[num_snd_max], U[num_snd_max], V[num_snd_max];
    FILE *in;

    if ( !(in = fopen(bg_fl, "r")) ) {
	fprintf(stderr, "Could not open background file %s.\n", bg_fl);
	exit(EXIT_FAILURE);
    }
    printf("Reading sounding file %s\n", bg_fl);
    for (n = 0; !feof(in) && n < num_snd_max; n++) {
	if ( fscanf(in, " %f %f %f", height + n, U + n, V + n) != 3 ) {
	    if ( ferror(in) ) {
		fprintf(stderr, "Could not read line %d of sounding file %s.\n",
			n, bg_fl);
		exit(EXIT_FAILURE);
	    }
	}
    }
    fclose(in);
    if ( n == num_snd_max ) {
	fprintf(stderr, "Too many sounding levels. "
		"Maximum allowed = %d\n", num_snd_max);
	exit(EXIT_FAILURE);
    }
    num_snd = n;

    /* Linearly interpolate sounding to analysis heights */
    USND[0] = U[0];
    VSND[0] = V[0];
    for (k = 1; k < num_z; k++) {
	USND[k] = VSND[k] = 0.0;
	z = z_min + k * z_inc;
	for (i = 0; i < num_snd; i++) {
	    if (height[i] > z) {
		double f;		/* Interpolation factor */

		f = (z - height[i - 1]) / (height[i] - height[i - 1]);
		USND[k] = U[i - 1] + (U[i] - U[i - 1]) * f;
		VSND[k] = V[i - 1] + (V[i] - V[i - 1]) * f;
		break;
	    }
	}
    }

    printf(" k        z     USND     VSND |  k        z     USND     VSND\n");
    for (k = 0; k < num_z / 2; k++) {
	z = z_min + k * z_inc;
	printf("%2d %8.1f %8.1f %8.1f | ", k, z, USND[k], VSND[k]);
	int k1 = k + num_z / 2;
	z = z_min + k1 * z_inc;
	printf("%2d %8.1f %8.1f %8.1f\n", k1, z, USND[k1], VSND[k1]);
    }
    if ( --k + num_z / 2 != num_z - 1 ) {
	printf("                              | %2d %8.1f %8.1f %8.1f\n",
		num_z - 1, z, USND[num_z - 1], VSND[num_z - 1]);
    }
    printf("\n");

    n = 0;
    for (i = 0; i < num_x; i++) {
	for (j = 0; j < num_y; j++) {
	    for (k = 0; k < num_z; k++) {
		n++;
		if (k >=  0 && k < num_z) {
		    ub[n] = USND[k];
		    vb[n] = VSND[k];
		}
	    }
	}
    }

}

/*
   Assign coverage values throughout grid.
   num_radar	- number of radars.
   num_ob_z	- number of vertical levels in observation grid.
   num_ob_y	- number of northward coordinates in observation grid.
   num_ob_x	- number of eastward coordinates in observation grid.
   z_ob		- array of vertical coordinates, meters, dimensioned [num_ob_z]
   y_ob		- array of northward coordinates, meters, dimensioned [num_ob_y]
   x_ob		- array of eastward coordinates, meters, dimensioned [num_ob_x]
   radz		- array of radar z coordinates, meters, dimensioned [num_radar]
   rady		- array of radar y coordinates, meters, dimensioned [num_radar]
   radx		- array of radar x coordinates, meters, dimensioned [num_radar]
   V_obs	- observed radial velocity, m/s,
   .		  dimensioned [num_radar * num_ob_z * num_ob_y * num_ob_x]
   refl_obs	- observed reflectivity,
   .		  dimensioned [num_radar * num_ob_z * num_ob_y * num_ob_x]
   V_obs_miss	- missing value for V_obs
   refl_obs_miss- missing value for reflectivity
   min_CBA	- minimum beam crossing angle for dual Doppler coverage,
   .		  degrees.

   If successful, global arrays coverage_bg, coverage_fil, SS, and max_dbz
   receive values.
 */

void init_wind3(int num_radar, size_t num_ob_z, size_t num_ob_y,
	size_t num_ob_x, float *z_ob, float *y_ob, float *x_ob, float *radz,
	double *rady, double *radx, float *V_obs, float *refl_obs,
	float V_obs_miss, float refl_obs_miss, double min_CBA)
{
    float ****vr;			/* Radial velocity field, V_obs */
    float ****refl;			/* Reflectivity, refl_obs */
    float ****ss;			/* SS field */
    int num_irng;			/* Number of radars with data at a
					   point */
    int irng[MAX_RADARS];		/* Array of indeces of radars with
					   data at a point.
					   For example, if num_irng = 2
					   and idcs = {0, 2, ...} then radars 0
					   and 2 are covering the point, radar 1
					   is not. */
    int count, cba_count;
    int n1, n2;
    int n, k, j, i;
    double s1, s2, s3, CBA;

    /*
       vr and refl hold V_Obs and refl_obs, respectively. The are both dimensioned
       [num_radar][num_ob_z][num_ob_y][num_ob_x].  Allocate and assign higher
       dimensions.
     */

    vr = calloc4f(num_radar, num_ob_z, num_ob_y, num_ob_x, V_obs, "vr");
    refl = calloc4f(num_radar, num_ob_z, num_ob_y, num_ob_x, refl_obs,
	    "reflectivity");
    ss = calloc4f(num_radar, num_ob_z, num_ob_y, num_ob_x, SS, "ss");

    printf("Computing coverage: ");
    count = 0;
    for (k = 0; k < num_ob_z; k++) {
	for (j = 0; j < num_ob_y; j++) {
	    for (i = 0; i < num_ob_x; i++) {
		maxdbz[i][j][k] = -INFINITY;

		/* Identify radars with velocity data at the point */
		for (n = 0; n < num_radar; n++) {
		    irng[n] = 0;
		}
		for (n = num_irng = 0; n < num_radar; n++) {
		    if ( vr[n][k][j][i] != V_obs_miss ) {
			irng[num_irng] = n;
			if ( refl[n][k][j][i] != refl_obs_miss
				&& maxdbz[i][j][k] < refl[n][k][j][i] ) {
			    maxdbz[i][j][k] = refl[n][k][j][i];
			}
			num_irng++;
		    }
		}

		/* Assign coverage flags to the point */
		if ( num_irng == 1 ){
		    /* Single Doppler */
		    n1 = irng[0];
		    ss[n1][k][j][i] = 0.0;
		    coverage_bg[i][j][k] = Get_Cvg_BG(n1);
		    coverage_fil[i][j][k] =  Get_Cvg_Fil(n1);
		} else if ( num_irng == 2 ){
		    /*
		       Two radars with data within range.
		       Check beam crossing angle.
		     */

		    n1 = irng[0];
		    n2 = irng[1];
		    s1 = dist2(x_ob[i], y_ob[j], radx[n1], rady[n1]);
		    s2 = dist2(x_ob[i], y_ob[j], radx[n2], rady[n2]);
		    s3 = dist2(radx[n2], rady[n2], radx[n1], rady[n1]);
		    CBA = Deg(acos((s1 * s1 + s2 * s2 - s3 * s3)
				/ (2 * s1 * s2)));
		    if ( fabs(CBA) >= min_CBA
			    && fabs(CBA) <= (180 - min_CBA) ) {
			/* Optimal beam crossing angle */
			ss[n1][k][j][i] = ss[n2][k][j][i] = 1.0;
			coverage_bg[i][j][k] = Get_Cvg_Opt_BG(2, irng);
			coverage_fil[i][j][k] = Get_Cvg_Opt_Fil(2, irng);
			count++;
		    } else {
			/* Sub-optimal beam crossing angle */
			ss[n1][k][j][i] = ss[n2][k][j][i] = 0.0;
			coverage_bg[i][j][k] = Get_Cvg_Sub_BG(2, irng);
			coverage_fil[i][j][k] = Get_Cvg_Sub_Fil(2, irng);
			count++;
		    }
		} else if ( num_irng == 3 ) {

		    /*
		       Three radars have data near the point.
		       Check beam crossing angle for all pairs.
		     */

		    for (cba_count = 0, n1 = 0; n1 < num_radar - 1; n1++) {
			for (n2 = n1 + 1; n2 < num_radar; n2++) {
			    s1 = dist2(x_ob[i], y_ob[j], radx[n1], rady[n1]);
			    s2 = dist2(x_ob[i], y_ob[j], radx[n2], rady[n2]);
			    s3 = dist2(radx[n2], rady[n2], radx[n1], rady[n1]);
			    CBA = Deg(acos((s1 * s1 + s2 * s2 - s3 * s3)
					/ (2 * s1 * s2)));
			    if ( fabs(CBA) >= min_CBA
				    || fabs(CBA) <= (180 - min_CBA) ) {
				cba_count++;
			    }
			}
		    }

		    if ( cba_count == 0 ) {
			/* No dual Dopper coverage */
			for (n1 = 0; n1 < num_radar; n1++) {
			    ss[n1][k][j][i] = 0.0;
			}
		    } else if ( cba_count < num_radar ) {
			/* Dual Doppler coverage, although not all radars */
			for (n1 = 0; n1 < num_radar; n1++) {
			    ss[n1][k][j][i] = 0.0;
			}
			coverage_bg[i][j][k] = 1;
			coverage_fil[i][j][k] = 1;
			count++;
		    } else if ( cba_count == num_radar ) {
			/* Dual Doppler coverage with all radars */
			for (n1 = 0; n1 < num_radar; n1++) {
			    ss[n1][k][j][i] = Get_SSeq_Trip(n1);
			}
			coverage_bg[i][j][k] = 1;
			coverage_fil[i][j][k] = 1;
			count++;
		    }
		} 
	    } 
	}
    }
    printf("%d verification points out of %zd total points.\n",
	    count, num_ob_z * num_ob_y * num_ob_x);

    /* If maxdbz is unknown, assume it is known to be 0.0. For consistency. */
    for (i = 0; i < num_x; i++) {
	for (j = 0; j < num_y; j++) {
	    for (k = 0; k < num_z; k++) {
		if ( !isfinite(maxdbz[i][j][k]) ) {
		    maxdbz[i][j][k] = 0.0;
		}
	    }
	}
    }

    /* Done with the higher dimensions */
    free4f(vr);
    free4f(refl);
    free4f(ss);
}

void ReadPyART(int num_radar, char opaws_fl_nms[MAX_RADARS][PATH_LEN_MAX],
	const char *refl_nm, const char *vt_nm, double min_CBA)
{
    int ncid[num_radar];		/* NetCDF identifiers for OPAWS files */
    jmp_buf nc_err;			/* Jump buffer for NetCDF errors */
    size_t num_ob_x, num_ob_y, num_ob_z;/* OPAWS grid dimensions */
    float *z_ob, *y_ob, *x_ob;		/* OPAWS grid coordinates */
    float *z_tmp, *y_tmp, *x_tmp;	/* Temporarily hold OPAWS grid
					   coordinates for comparison */
    float z_anal, y_anal, x_anal;	/* Analysis grid coordinates */
    int i, j, k, n = -1;
    int o1, o2;				/* Observation index */
    size_t sz;				/* Size of a new allocation */
    struct GeogProj geog_proj;          /* Map projection */
    double d_lat;                       /* Distance in great circle degrees from
					   grid origin to projection parallel */
    int radcount[MAX_RADARS] = {0};
    int retval;
    float z_diff_max = 10.0;		/* Maximum allowed difference between
					   OPAWS and analysis z coordinate */
    float y_diff_max = 10.0;		/* Maximum allowed difference between
					   OPAWS and analysis y coordinate */
    float x_diff_max = 10.0;		/* Maximum allowed difference between
					   OPAWS and analysis x coordinate */
    int match;				/* If true, grids match */
    float *azims_miss, *elevs_miss,
	  *V_obs_miss, *refl_obs_miss;	/* Missing values */

    printf("Py-ART grids for analyis and all radars must match to within "
	    "(dx, dy, dz) < (%g %g %g) meters.\n",
	    x_diff_max, y_diff_max, z_diff_max);

    /*
       Geographic projection to compute grid coordinates of radars,
       which are absent from OPAWS NetCDF file. This is should be the
       same projection OPAWS uses, although there is no way this program
       can verify such.
     */

    d_lat = (y_max - y_min) / NMILE / 60.0 * RAD_DEG / 6.0;
    if ( !(GeogProjSetLambertConfConic(grid_lon, grid_lat,
		    grid_lat - d_lat, grid_lat + d_lat, &geog_proj)) ) {
	fprintf(stderr, "Failed to set map projection.\n");
	exit(EXIT_FAILURE);
    }

    /* Set a jump buffer. Come back here if a NetCDF error occurs. */
    if ( setjmp(nc_err) == NNCDF_ERROR ) {
	fprintf(stderr, "Could not read Py-ART file%s.\n",
		(n == -1) ? "s" : opaws_fl_nms[n]);
	exit(EXIT_FAILURE);
    }

    /* Open the OPAWS files. Ensure they contain data for one time only. */
    for (n = 0; n < num_radar; n++) {
	ncid[n] = NNC_Open(opaws_fl_nms[n], nc_err);
	if ( NNC_Inq_Dim(ncid[n], "time", nc_err) != 1 ) {
	    fprintf(stderr, "%s has data for more than one time.\n",
		    opaws_fl_nms[n]);
	    exit(EXIT_FAILURE);
	}
    }

    /* Get specifications for observation grid from first radar */
    z_ob = y_ob = x_ob = NULL;
    n = 0;
    num_ob_z = NNC_Inq_Dim(ncid[n], "z", nc_err);
    num_ob_y = NNC_Inq_Dim(ncid[n], "y", nc_err);
    num_ob_x = NNC_Inq_Dim(ncid[n], "x", nc_err);
    z_ob = falloc(num_ob_z, "z_ob");
    y_ob = falloc(num_ob_y, "y_ob");
    x_ob = falloc(num_ob_x, "x_ob");
    NNC_Get_Var_Float(ncid[n], "z", z_ob, nc_err);
    NNC_Get_Var_Float(ncid[n], "y", y_ob, nc_err);
    NNC_Get_Var_Float(ncid[n], "x", x_ob, nc_err);

    /* Convert m to m */
    for (k = 0; k < num_ob_z; k++) {
	z_ob[k] *= 1.0;
    }
    for (j = 0; j < num_ob_y; j++) {
	y_ob[j] *= 1.0;
    }
    for (i = 0; i < num_ob_x; i++) {
	x_ob[i] *= 1.0;
    }

    /*
       Make sure observation grid is compatible with analysis grid.
       Note that the analysis grid could have additional points above,
       to the east, or to the north.
     */

    match = 1;
    for (k = 0; k < num_ob_z; k++) {
	z_anal = z_min + k * z_inc;
	if ( fabs(z_ob[k] - z_anal) > z_diff_max ) {
	    fprintf(stderr, "Z grid mismatch: z_ob[%d] = %g. "
		    "z_analz = %g. Diff = %g.\n",
		    k, z_ob[k], z_anal, fabs(z_ob[k] - z_anal));
	    match = 0;
	}
    }
    for (j = 0; j < num_ob_y; j++) {
	y_anal = y_min + j * y_inc;
	if ( fabs(y_ob[j] - y_anal) > y_diff_max ) {
	    fprintf(stderr, "Y grid mismatch: y_ob[%d] = %g. "
		    "y_analz = %g. Diff = %g.\n",
		    j, y_ob[j], y_anal, fabs(y_ob[j] - y_anal));
	    match = 0;
	}
    }

    for (i = 0; i < num_ob_x; i++) {
	x_anal = x_min + i * x_inc;
	if ( fabs(x_ob[i] - x_anal) > x_diff_max ) {
	    fprintf(stderr, "X grid mismatch: x_ob[%d] = %g. "
		    "x_analz = %g. Diff = %g.\n",
		    i, x_ob[i], x_anal, fabs(x_ob[i] - x_anal));
	    match = 0;
	}
    }
    if ( !match ) {
	fprintf(stderr, "Grid in %s does not match analysis.\n",
		opaws_fl_nms[n]);
	exit(EXIT_FAILURE);
    }

    /* Make sure grids for other radars match the first */
    z_tmp = falloc(num_ob_z, "z_ob");
    y_tmp = falloc(num_ob_y, "y_ob");
    x_tmp = falloc(num_ob_x, "x_ob");
    for (n = 1; n < num_radar; n++) {
	match = 1;
	if ( NNC_Inq_Dim(ncid[n], "z", nc_err) != num_ob_z ) {
	    fprintf(stderr, "Number of z grid points in %s does not "
		    "match observation grid. Should have been %zd. "
		    "Got %zd\n", opaws_fl_nms[n], num_ob_z, 
		    NNC_Inq_Dim(ncid[n], "z", nc_err));
	    exit(EXIT_FAILURE);
	}
	if ( NNC_Inq_Dim(ncid[n], "y", nc_err) != num_ob_y ) {
	    fprintf(stderr, "Number of y grid points in %s does not "
		    "match observation grid. Should have been %zd. "
		    "Got %zd\n", opaws_fl_nms[n], num_ob_y, 
		    NNC_Inq_Dim(ncid[n], "y", nc_err));
	    exit(EXIT_FAILURE);
	}
	if ( NNC_Inq_Dim(ncid[n], "x", nc_err) != num_ob_x ) {
	    fprintf(stderr, "Number of x grid points in %s does not "
		    "match observation grid. Should have been %zd. "
		    "Got %zd\n", opaws_fl_nms[n], num_ob_x, 
		    NNC_Inq_Dim(ncid[n], "x", nc_err));
	    exit(EXIT_FAILURE);
	}
	NNC_Get_Var_Float(ncid[n], "z", z_tmp, nc_err);
	NNC_Get_Var_Float(ncid[n], "y", y_tmp, nc_err);
	NNC_Get_Var_Float(ncid[n], "x", x_tmp, nc_err);

	/* Convert m to m */
	for (k = 0; k < num_ob_z; k++) {
	    z_tmp[k] *= 1.0;
	}
	for (j = 0; j < num_ob_y; j++) {
	    y_tmp[j] *= 1.0;
	}
	for (i = 0; i < num_ob_x; i++) {
	    x_tmp[i] *= 1.0;
	}

	for (k = 0; k < num_ob_z; k++) {
	    if ( fabs(z_ob[k] - z_tmp[k]) > z_diff_max ) {
		fprintf(stderr, "%s z grid off. z[%d] = %g instead of %g.\n",
			opaws_fl_nms[n], k, z_tmp[k], z_ob[k]);
		match = 0;
	    }
	}
	for (j = 0; j < num_ob_y; j++) {
	    if ( fabs(y_ob[j] - y_tmp[j]) > y_diff_max ) {
		fprintf(stderr, "%s y grid off. y[%d] = %g instead of %g.\n",
			opaws_fl_nms[n], j, y_tmp[j], y_ob[j]);
		match = 0;
	    }
	}
	for (i = 0; i < num_ob_x; i++) {
	    if ( fabs(x_ob[i] - x_tmp[i]) > x_diff_max ) {
		fprintf(stderr, "%s x grid off. x[%d] = %g instead of %g.\n",
			opaws_fl_nms[n], i, x_tmp[i], x_ob[i]);
		match = 0;
	    }
	}
	if ( !match ) {
	    fprintf(stderr, "Observation grid in %s does not match.\n",
		    opaws_fl_nms[n]);
	    exit(EXIT_FAILURE);
	}
    }
    free(z_tmp);
    free(y_tmp);
    free(x_tmp);

    /*
       Initialize arrays of observed values to accomodate all grid
       points for all radars. MM and associated allocations might be
       reduced later as missing or unusable values are discarded.
     */

    MM = num_radar * num_ob_z * num_ob_y * num_ob_x;
    zpos = falloc(MM, "zpos");
    ypos = falloc(MM, "ypos");
    xpos = falloc(MM, "xpos");
    map = ialloc(MM, "map");
    azims = falloc(MM, "azims");
    elevs = falloc(MM, "elevs");
    V_obs = falloc(MM, "V_obs");
    refl_obs = falloc(MM, "refl_obs");
    term_vel = falloc(MM, "term_vel");
    SS = falloc(MM, "SS");


    /*
       Assign observation grid information to global arrays.
       Note that global nn array is dimensioned [num_x][num_y][num_z],
       which are the dimensions of the analysis grid.  Before looping
       through observation grid, make sure assignment to map array below
       stays in bounds.
     */

    if ( num_ob_z > num_z ) {
	fprintf(stderr, "Observation domain cannot have more z levels (%zd) "
		"than analysis domain (%d)\n", num_ob_z, num_z);
	exit(EXIT_FAILURE);
    }
    if ( num_ob_y > num_y ) {
	fprintf(stderr, "Observation domain cannot have more y levels (%zd) "
		"than analysis domain (%d)\n", num_ob_y, num_y);
	exit(EXIT_FAILURE);
    }
    if ( num_ob_y > num_y ) {
	fprintf(stderr, "Observation domain cannot have more y levels (%zd) "
		"than analysis domain (%d)\n", num_ob_y, num_y);
	exit(EXIT_FAILURE);
    }
    o1 = 0;
    for (i = 0; i < num_ob_x; i++) {
	for (j = 0; j < num_ob_y; j++) {
	    for (k = 0; k < num_ob_z; k++) {
		map[o1] = nn[i][j][k];
		o1++;
	    }
	}
    }

    /* Read data */
    printf("Reading Py-ART file");
    for (n = 0; n < num_radar; n++) {
	float radlon, radlat; 		/* Radar longitude, latitude */
	int offset;			/* Offset to start of data for
					   current radar in a data array. */

	printf("%s\n", opaws_fl_nms[n]);
	NNC_Get_Var_Float(ncid[n], "radar_longitude", &radlon, nc_err);
	NNC_Get_Var_Float(ncid[n], "radar_latitude", &radlat, nc_err);
	NNC_Get_Var_Float(ncid[n], "radar_altitude", radz + n, nc_err);
	radz[n] *= 1000.0;
	if ( !GeogProjLonLatToXY(radlon * RAD_DEG, radlat * RAD_DEG,
		    radx + n, rady + n, &geog_proj) ) {
	    fprintf(stderr, "Could not compute map coordinates of radar %d "
		    "in current projection.\n", n);
	    exit(EXIT_FAILURE);
	}

	offset = n * num_ob_z * num_ob_y * num_ob_x;
	NNC_Get_Var_Float(ncid[n], "AZ", azims + offset, nc_err);
	NNC_Get_Var_Float(ncid[n], "EL", elevs + offset, nc_err);
	NNC_Get_Var_Float(ncid[n], vt_nm, V_obs + offset, nc_err);
	NNC_Get_Var_Float(ncid[n], refl_nm, refl_obs + offset, nc_err);
    }
    printf("\n");
    printf("Radar positions: ");
    for (n = 0; n < num_radar; n++) {
	printf(" (%g,%g)", radx[n], rady[n]);
    }
    printf("\n");

    /*
       Get flags for missing values for the first radar.  Assume, perhaps
       naively, that all radars represent missing values the same way.
     */

    n = 0;
    azims_miss = elevs_miss = V_obs_miss = refl_obs_miss = NULL;
    azims_miss = NNC_Get_Att_Float(ncid[n], "AZ", "missing_value", nc_err);
    elevs_miss = NNC_Get_Att_Float(ncid[n], "EL", "missing_value", nc_err);
    V_obs_miss = NNC_Get_Att_Float(ncid[n], vt_nm, "missing_value", nc_err);
    refl_obs_miss = NNC_Get_Att_Float(ncid[n], refl_nm, "missing_value",
	    nc_err);

    init_wind3(num_radar, num_ob_z, num_ob_y, num_ob_x, z_ob, y_ob, x_ob,
	    radz, rady, radx, V_obs, refl_obs, V_obs_miss[0], refl_obs_miss[0],
	    min_CBA);

    /*
       Discard undesirable points. o1 is the index of the last desired point.
       o2 is the index of the current point in the loop. If observation at
       o2 is good, move it to o1, and increment o1.
     */

    for (n = o1 = o2 = 0; n < num_radar; n++) {
	for (k = 0; k < num_ob_z; k++) {
	    obs_index[k][n] = o1;
	    for (j = 0; j < num_ob_y; j++) {
		for (i = 0; i < num_ob_x; i++) {
		    if ( fabs(V_obs[o2]) < 200 ) {
			radcount[n]++;
			zpos[o1] = z_ob[k];
			ypos[o1] = y_ob[j];
			xpos[o1] = x_ob[i];
			map[o1] = nn[i][j][k];
			azims[o1] = azims[o2];
			elevs[o1] = elevs[o2];
			V_obs[o1] = V_obs[o2];
			refl_obs[o1] = refl_obs[o2];
			term_vel[o1] = CalcTermVel(refl_obs[o1], z_ob[k]);
			SS[o1] = SS[o2];
			o1++;
		    }
		    o2++;
		}
	    }
	}
    }
    MM = o1;

    FREE(azims_miss);
    FREE(elevs_miss);
    FREE(V_obs_miss);
    FREE(refl_obs_miss);

    for (n = 0; n < num_radar; n++) {
	retval = nc_close(ncid[n]);
	if ( retval != NC_NOERR ) {
	    fprintf(stderr, "Could not close %s.\n%s\n",
		    opaws_fl_nms[n], nc_strerror(retval));
	    exit(EXIT_FAILURE);
	}
    }
    free(z_ob);
    free(y_ob);
    free(x_ob);

    /* Adjust allocations to observation count. */
    sz = MM * sizeof(float);
    if ( !(V_obs = (float *)realloc(V_obs, sz)) ) {
	fprintf(stderr, "Could not reallocate V_obs\n");
	exit(EXIT_FAILURE);
    }
    if ( !(xpos = (float *)realloc(xpos, sz)) ) {
	fprintf(stderr, "Could not reallocate xpos\n");
	exit(EXIT_FAILURE);
    }
    if ( !(ypos = (float *)realloc(ypos, sz)) ) {
	fprintf(stderr, "Could not reallocate ypos\n");
	exit(EXIT_FAILURE);
    }
    if ( !(zpos = (float *)realloc(zpos, sz)) ) {
	fprintf(stderr, "Could not reallocate zpos\n");
	exit(EXIT_FAILURE);
    }
    if ( !(elevs = (float *)realloc(elevs, sz)) ) {
	fprintf(stderr, "Could not reallocate elevs\n");
	exit(EXIT_FAILURE);
    }
    if ( !(azims = (float *)realloc(azims, sz)) ) {
	fprintf(stderr, "Could not reallocate azims\n");
	exit(EXIT_FAILURE);
    }
    if ( !(refl_obs = (float *)realloc(refl_obs, sz)) ) {
	fprintf(stderr, "Could not reallocate "
		"refl_obs\n");
	exit(EXIT_FAILURE);
    }
    if ( !(term_vel = (float *)realloc(term_vel, sz)) ) {
	fprintf(stderr, "Could not reallocate "
		"term_vel\n");
	exit(EXIT_FAILURE);
    }
    sz = MM * sizeof(int);
    if ( !(map = (int *)realloc(map, sz)) ) {
	fprintf(stderr, "Could not reallocate map\n");
	exit(EXIT_FAILURE);
    }
    for (n = 0; n < num_radar; n++) {
	printf("radar%d %d obs. ", n, radcount[n]);
    }
    printf("%d total obs.\n", MM);
}

double CalcTermVel(double refl, double z) {

    // Method provided by Mike Biggerstaff and Dan Betten

    double A,B,fallspeed,rho;
    double frz=4500;

    if (z < frz) {

	if (refl < 55)      { A=-2.6;  B=.0107; } // rain
	else if (refl < 60) { A=-2.5;  B=.013;  } // graupel
	else                { A=-3.95; B=.0148; } // hail

    }

    else {

	if (refl < 33)       { A=-.817; B=.0063; } // ice
	else if (refl < 49)  { A=-2.5;  B=.013; } // graupel
	else                 { A=-3.95; B=.0148; } // hail

    }

    rho=exp(-z/10000);
    fallspeed=A*pow(10, refl*B)*pow(1.2/rho,0.4);

    return fallspeed;

}

// Continuity cost function
double CalcDiv (double xp[]) {
    double DIV=0;
    double z;
    int n=0, count=0;
    double anel_term;

    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {

		n=nn[i][j][k];

		if ( i > 0 && i < num_dx ) {
		    dudx[i][j][k]=(xp[n+num_z*num_y]-xp[n-num_z*num_y])/(2*x_inc);
		    dvdx[i][j][k]=(xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc);
		    dwdx[i][j][k]=(xp[n+num_z*num_y+2*N]-xp[n-num_z*num_y+2*N])/(2*x_inc);
		}
		if ( j > 0 && j < num_dy ) {
		    dudy[i][j][k]=(xp[n+num_z]-xp[n-num_z])/(2*y_inc);
		    dvdy[i][j][k]=(xp[n+num_z+N]-xp[n-num_z+N])/(2*y_inc);
		    dwdy[i][j][k]=(xp[n+num_z+2*N]-xp[n-num_z+2*N])/(2*y_inc);
		}
		if ( k > 0 && k < num_dz ) {
		    dudz[i][j][k]=(xp[n+1]-xp[n-1])/(2*z_inc);
		    dvdz[i][j][k]=(xp[n+1+N]-xp[n-1+N])/(2*z_inc);
		    dwdz[i][j][k]=(xp[n+1+2*N]-xp[n-1+2*N])/(2*z_inc);
		}
		if (k<num_dz) {
		    dwdz3[i][j][k]=(xp[n+1+2*N]-xp[n+2*N])/z_inc;
		}

		if ( i > 0 && j > 0 && i < num_dx && j < num_dy )
		    vort[i][j][k]=(xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc)
			-(xp[n+num_z]-xp[n-num_z])/(2*y_inc);

	    }
	}
    }
    n=0;

    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {

		n=nn[i][j][k];

		if ( i > 0 && i < num_dx
			&& j > 0 && j < num_dy
			&& k > -1 && k < num_dz ) {
		    z=z_min+k*z_inc;

		    // check if in buffer domain
		    if (z>=z_min) {
			if ( i < 0 || i > num_dx
				|| j < 0 || j > num_dy
				|| k < 0 || k > num_dz) {
			    C2=C2d;
			} else {
			    C2=C2a;
			}

			// if anelastic approx. used, compute for extra term
			if ( anel==1 ) {
			    anel_term = 1.0 / ((rho[k]+rho[k+1])/2)
				* (rho[k+1]-rho[k])/(z_inc);
			} 
			else {
			    anel_term=0;
			}
			DIV += C2 * pow((dudx[i][j][k]+dudx[i][j][k+1]) / 2
				+ (dvdy[i][j][k]+dvdy[i][j][k+1]) / 2
				+ beta * (dwdz3[i][j][k]+xp[n+2*N]*anel_term)
				,2);
			div2[i][j][k] = (dudx[i][j][k]+dudx[i][j][k+1]) / 2
			    + (dvdy[i][j][k]+dvdy[i][j][k+1]) / 2
			    + beta * (dwdz3[i][j][k]+xp[n+2*N] * anel_term);
			count++;
		    }
		}
	    }
	}
    }
    return DIV;
}

// Vorticity equation cost function
double CalcVortTend(double xp[]) {
    double z,vort_tend=0,dzetadt,dzetadx,dzetady,dzetadz;
    int n=0;
    int count=0;
    double tilting, stretching;

    n=0;
    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {
		n++;
		if (j==(num_dy-1)&&i>1&&i<(num_dx-1)&&k>0&&k<num_dz) {
		    dzetadx=(vort[i+1][j][k]-vort[i-1][j][k])/(2*x_inc);
		    dzetady=(dvdx[i][j+1][k]-dvdx[i][j-1][k])/(2*y_inc)-(xp[n+num_z]-2*xp[n]+xp[n-num_z])/(y_inc*y_inc);
		    dzetadz=(vort[i][j][k+1]-vort[i][j][k-1])/(2*z_inc);
		    tilting=dvdz[i][j][k]*dwdx[i][j][k]-dudz[i][j][k]*dwdy[i][j][k];
		    stretching=-(vort[i][j][k])*dwdz[i][j][k];

		}

		if (i>1&&i<(num_dx-1)&&j>1&&j<(num_dy-1)&&k>0&&k<num_dz) {

		    z=z_min+k*z_inc;

		    // local derivative accounts only for contributions from advection (assumes steady-state wind field, unlike in Potvin et al. 2012)

		    dzetadt=VT*((dudy[i][j+1][k]-dudy[i][j-1][k])/(2*y_inc)-(dvdx[i][j+1][k]-dvdx[i][j-1][k])/(2*y_inc)) +
			UT*((dudx[i][j+1][k]-dudx[i][j-1][k])/(2*y_inc)-(dvdx[i+1][j][k]-dvdx[i-1][j][k])/(2*x_inc));

		    dzetadx=(vort[i+1][j][k]-vort[i-1][j][k])/(2*x_inc);
		    dzetady=(vort[i][j+1][k]-vort[i][j-1][k])/(2*y_inc);
		    dzetadz=(vort[i][j][k+1]-vort[i][j][k-1])/(2*z_inc);
		    tilting=dvdz[i][j][k]*dwdx[i][j][k]-dudz[i][j][k]*dwdy[i][j][k];
		    //stretching=-(vort[i][j][k])*dwdz[i][j][k]; // Boussinesq
		    stretching=(vort[i][j][k])*(dudx[i][j][k]+dvdy[i][j][k]); // anelastic
		    DzetaDt[i][j][k]=dzetadt+xp[n]*dzetadx+xp[n+N]*dzetady+xp[n+2*N]*dzetadz+tilting+stretching;

		    if (z >= z_min + z_inc) {
			if ( i < 0 || i > num_dx
				|| j < 0 || j > num_dy
				|| k < 0 || k > num_dz )
			    C3=C3d;
			else
			    C3=C3a;
			vort_tend+=C3*DzetaDt[i][j][k]*DzetaDt[i][j][k];
		    }

		    count++;
		}
	    }
	}
    }
    return vort_tend;
}

// Smoothness cost function
double CalcSmooth (double xp[]) {
    int n=0;
    double smoothness=0,z;

    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {

		n++;

		if (i>0&&i<=(num_dx-0)&&j>0&&j<=(num_dy-0)&&k>0&&k<=(num_dz-0)) {
		    if (laplace==1 && i>0&&i<num_dx&&j>0&&j<num_dy&&k>0&&k<num_dz) {

			dudx2[i][j][k]=(xp[n+num_z*num_y]-2*xp[n]+xp[n-num_z*num_y])/(x_inc*x_inc);
			dvdx2[i][j][k]=(xp[n+N+num_z*num_y]-2*xp[n+N]+xp[n+N-num_z*num_y])/(x_inc*x_inc);
			dwdx2[i][j][k]=(xp[n+2*N+num_z*num_y]-2*xp[n+2*N]+xp[n+2*N-num_z*num_y])/(x_inc*x_inc);

			dudy2[i][j][k]=(xp[n+num_z]-2*xp[n]+xp[n-num_z])/(y_inc*y_inc);
			dvdy2[i][j][k]=(xp[n+N+num_z]-2*xp[n+N]+xp[n+N-num_z])/(y_inc*y_inc);
			dwdy2[i][j][k]=(xp[n+2*N+num_z]-2*xp[n+2*N]+xp[n+2*N-num_z])/(y_inc*y_inc);

			dudz2[i][j][k]=(xp[n+1]-2*xp[n]+xp[n-1])/(z_inc*z_inc);
			dvdz2[i][j][k]=(xp[n+N+1]-2*xp[n+N]+xp[n+N-1])/(z_inc*z_inc);
			dwdz2[i][j][k]=(xp[n+2*N+1]-2*xp[n+2*N]+xp[n+2*N-1])/(z_inc*z_inc);

			z=z_min+k*z_inc;

			if (z<cutoff) sss=ss; else sss=1.0;

			if (z >= z_min) {

			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C4=C4d; else C4=C4a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C5=C5d; else C5=C5a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C6=C6d; else C6=C6a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C7=C7d; else C7=C7a;

			    smoothness+=C5*sss*(dudz2[i][j][k]*dudz2[i][j][k]+dvdz2[i][j][k]*dvdz2[i][j][k]);
			    smoothness+=C6*sss*(dwdz2[i][j][k]*dwdz2[i][j][k]);
			    smoothness+=C7*sss*(dwdx2[i][j][k]*dwdx2[i][j][k]+dwdy2[i][j][k]*dwdy2[i][j][k]);
			    smoothness+=C4*sss*(dudx2[i][j][k]*dudx2[i][j][k]+dudy2[i][j][k]*dudy2[i][j][k]+dvdx2[i][j][k]*dvdx2[i][j][k]+dvdy2[i][j][k]*dvdy2[i][j][k]);

			}

		    } // laplace=1

		    else if (laplace==0 && i>0 && j>0 && k>0) {

			dwdx2[i][j][k]=(xp[n+2*N]-xp[n-num_z*num_y+2*N])/x_inc;
			dwdy2[i][j][k]=(xp[n+2*N]-xp[n-num_z+2*N])/y_inc;
			dwdz2[i][j][k]=(xp[n+2*N]-xp[n-1+2*N])/z_inc;

			dudx2[i][j][k]=(xp[n]-xp[n-num_z*num_y])/x_inc;
			dudy2[i][j][k]=(xp[n]-xp[n-num_z])/y_inc;

			dvdx2[i][j][k]=(xp[n+N]-xp[n-num_z*num_y+N])/x_inc;
			dvdy2[i][j][k]=(xp[n+N]-xp[n-num_z+N])/y_inc;

			dudz2[i][j][k]=(xp[n]-xp[n-1])/z_inc;
			dvdz2[i][j][k]=(xp[n+N]-xp[n-1+N])/z_inc;

			z=z_min+k*z_inc;

			if (z<=cutoff) sss=ss; else sss=1.0;

			if (z >= z_min) {

			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C4=C4d; else C4=C4a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C5=C5d; else C5=C5a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C6=C6d; else C6=C6a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz) C7=C7d; else C7=C7a;

			    if (z<=cutoff) {
				smoothness+=C5*sss*(dudz2[i][j][k]*dudz2[i][j][k]+dvdz2[i][j][k]*dvdz2[i][j][k]);
				smoothness+=C6*sss*(dwdz2[i][j][k]*dwdz2[i][j][k]);
			    }
			    else {
				smoothness+=C5*sss*(dudz2[i][j][k]*dudz2[i][j][k]+dvdz2[i][j][k]*dvdz2[i][j][k]);
				smoothness+=C6*sss*(dwdz2[i][j][k]*dwdz2[i][j][k]);
			    }
			    smoothness+=C7*sss*(dwdx2[i][j][k]*dwdx2[i][j][k]+dwdy2[i][j][k]*dwdy2[i][j][k]);
			    smoothness+=C4*sss*(dudx2[i][j][k]*dudx2[i][j][k]+dudy2[i][j][k]*dudy2[i][j][k]+dvdx2[i][j][k]*dvdx2[i][j][k]+dvdy2[i][j][k]*dvdy2[i][j][k]);

			}
		    } // laplace=0
		}

	    }
	}
    }
    return smoothness;
}

// Background cost function
double CalcBG (double xp[]) {
    int n=0;
    double BG=0;

    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {
		n++;
		if ( i >= 2 && i < num_x-2
			&& j >= 2 && j < num_y-2
			&& k >= 0 && k < num_z ) {
		    if ( coverage_bg[i][j][k] == 0) {
			if (maxdbz[i][j][k]>10  && k<=weak_height){
			    BG+=C8a*(pow(xp[n]-ub[n],2)+pow(xp[n+N]-vb[n],2))/100.0;}
			else {
			    BG+=C8a*(pow(xp[n]-ub[n],2)+pow(xp[n+N]-vb[n],2));}
		    }
		    if ( (i <= 1 || i >= num_x-1
				|| j <= 1 || j >= num_y-1)
			    && (k >= 0 && k < num_z )) {
			if ( coverage_bg[i][j][k] == 0) {
			    BG+=C8a*(pow(xp[n]-ub[n],2)+pow(xp[n+N]-vb[n],2));}
		    }

		}

	    }
	}
    }
    return BG;
}

// Compute J (cost function)
double CalcCost (double xp[]) {
    double V_mod;			/* Model Vr */
    double uu, vv, ww;
    int i, n;
    double toterror, toterrorz, totwindz, zcount;
    double totwind;
    double J;				/* Return value = cost */
    double *JJ, *totwind2, *toterrorr;
    double meanwind, bgtotal;
    double elev_cutoff = 60.0;		/* Maximum tilt used in the obs error
					   calculation */

    JJ = dalloc(MM, "JJ");
    totwind2 = dalloc(MM, "totwind2");
    toterrorr = dalloc(MM, "toterrorr");

    zcount = toterrorz = totwindz = 0.0;
    obstotal=masstotal=vorttotal=bgtotal= 0.0;
    maxobserror=0.0;
    toterror = 0.0;
    J = 0.0;
    Jtotal=0.0;

    if (C2b>0)
	masstotal= CalcDiv(xp);		
    J += masstotal;		// Continuity cost function
    if (C3b>0)
	vorttotal = CalcVortTend(xp);
    J += vorttotal;		// Vorticity cost function
    if (C4b>0 || C5b>0 || C6b>0 || C7b>0)
    smoothtotal = CalcSmooth(xp);
    J += smoothtotal;		// Smoothness cost function
    if (C8b>0)
	bgtotal = CalcBG(xp);		// Background cost function
    J += bgtotal;
    totwind = 0.0;
    for (i = 0; i < MM; i++) {		// Data cost function
	totwind2[i] = V_obs[i] * V_obs[i];
	n = map[i];
	uu = xp[n];
	vv = xp[n+N];
	ww = xp[n+2*N];
	V_mod = uu * cos(Rad(elevs[i])) * sin(Rad(azims[i]))
	    + vv * cos(Rad(elevs[i])) * cos(Rad(azims[i]))
	    + (ww + term_vel[i]) * sin(Rad(elevs[i]));
	if (Vr_error==1)
	    vr_retr[i] = V_mod;
	JJ[i] = W1a[i] * ((V_mod - V_obs[i]) * (V_mod - V_obs[i]));
	toterrorr[i] = SS[i] * (V_mod - V_obs[i]) * (V_mod - V_obs[i]);

	// Store for use in CalcGrad:
	p[i]=2*(V_mod-V_obs[i])*cos(Rad(elevs[i]))*sin(Rad(azims[i]));
	p[i+MM]=2*(V_mod-V_obs[i])*cos(Rad(elevs[i]))*cos(Rad(azims[i]));
	p[i+2*MM]=2*(V_mod-V_obs[i])*sin(Rad(elevs[i]));
	obstotal += JJ[i];
    J += JJ[i];
	totwind += totwind2[i];
	toterror += toterrorr[i];
	// Estimate obs error for elevation below 5 degrees
	if ( elevs[i] <= elev_cutoff ) {
	    totwindz += totwind2[i];
	    toterrorz += toterrorr[i];
	    if (maxobserror<toterrorr[i]) maxobserror=toterrorr[i];
	    zcount +=1;
	}
    }
    rmsVr = meanwind = sqrt(totwind / MM);
    meanobserror = sqrt(toterrorz/totwindz);
    Jtotal=J-bgtotal;
    if (docount==1)
	printf("Analyzed obs = %d. Mean Vr = %g\n", MM, meanwind);
    docount=0;
    free(JJ);
    free(totwind2);
    free(toterrorr);
    return J;
}

// Adjoint
void CalcGrad(double* xp, double* grad)
{
    double z;
    int n=0;
    double anel_term;
    int ii, jj, i, j, k;

    for (ii=1; ii <= 3 * N; ii++) {
	grad[ii]=0;
    }

    // *** Adjoint for data cost function

    for (jj = 0; jj < MM; jj++) {
	n = map[jj];
	grad[n] += W1a[jj] * p[jj];
	grad[n+N] += W1a[jj] * p[jj + MM];
	grad[n+2*N] += W1a[jj] * p[jj +2*MM];
    }

    // *** Adjoint for CalcDiv

    if (C2b>0) {
	n=0;
	for (i=0; i < num_x; i++) {
	    for (j=0; j < num_y; j++) {
		for (k=0; k < num_z; k++) {
		    n++;
		    if (i>0&&i<num_dx&&j>0&&j<num_dy&&k>-1&&k<num_dz) {
			z=z_min+k*z_inc;
			if (z>=z_min) {
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C2=C2d;
			    else
				C2=C2a;
			    if (anel==1)
				anel_term=1.0/((rho[k]+rho[k+1])/2)*(rho[k+1]-rho[k])/z_inc;
			    else
				anel_term=0;
			    grad[n+num_z*num_y]+=2*C2*(div2[i][j][k]/(4*x_inc));
			    grad[n-num_z*num_y]-=2*C2*(div2[i][j][k]/(4*x_inc));
			    grad[n+num_z+N]+=2*C2*(div2[i][j][k]/(4*y_inc));
			    grad[n-num_z+N]-=2*C2*(div2[i][j][k]/(4*y_inc));

			    grad[n+1+num_z*num_y]+=2*C2*(div2[i][j][k]/(4*x_inc));
			    grad[n+1-num_z*num_y]-=2*C2*(div2[i][j][k]/(4*x_inc));
			    grad[n+1+num_z+N]+=2*C2*(div2[i][j][k]/(4*y_inc));
			    grad[n+1-num_z+N]-=2*C2*(div2[i][j][k]/(4*y_inc));

			    grad[n+1+2*N]+=2*C2*beta*div2[i][j][k]/z_inc;
			    grad[n+2*N]+=2*C2*beta*div2[i][j][k]*(anel_term-1/z_inc);

			}
		    }
		}
	    }
	}
    }

    // *** Adjoint for CalcVortTend

    if (C3b>0) {
	n=0;
	for (i=0; i < num_x; i++) {
	    for (j=0; j < num_y; j++) {
		for (k=0; k < num_z; k++) {
		    n++;
		    if (i>1&&i<(num_dx-1)&&j>1&&j<(num_dy-1)&&k>0&&k<num_dz) {
			z=z_min+k*z_inc;
			if (z >= z_min + z_inc) {
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C3=C3d;
			    else
				C3=C3a;
			    grad[n]+=-VT/(2*y_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+N]+=UT/(2*x_inc*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+num_z]+=UT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+num_z+N]+=-VT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y-num_z]+=UT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y-num_z+N]+=-VT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+num_z]+=-UT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+num_z+N]+=VT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y-num_z]+=-UT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y-num_z+N]+=VT/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+2*num_z]+=VT/(4*y_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-2*num_z]+=VT/(4*y_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+2*num_z*num_y+N]+=-UT/(4*x_inc*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-2*num_z*num_y+N]+=-UT/(4*x_inc*x_inc)*2*DzetaDt[i][j][k]*C3;

			    //horizontal advection

			    grad[n]+=((xp[n+2*num_z*num_y+N]-2*xp[n+N]+xp[n-2*num_z*num_y+N])/(4*x_inc*x_inc)-
				    ((xp[n+num_z*num_y+num_z]-xp[n+num_z*num_y-num_z])-(xp[n-num_z*num_y+num_z]-
					xp[n-num_z*num_y-num_z]))/(4*x_inc*y_inc)+xp[n+N]/(2*y_inc*y_inc))*2*DzetaDt[i][j][k]*C3;
			    grad[n+N]+=(-(xp[n+2*num_z]-2*xp[n]+xp[n-2*num_z])/(4*y_inc*y_inc)+
				    (xp[n+num_z*num_y+num_z+N]-xp[n-num_z*num_y+num_z+N]-(xp[n+num_z*num_y-num_z+N]-
											  xp[n-num_z*num_y-num_z+N]))/(4*x_inc*y_inc)-xp[n]/(2*x_inc*x_inc))*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+num_z]+=-xp[n]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+num_z+N]+=xp[n+N]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y-num_z]+=-xp[n]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y-num_z+N]+=xp[n+N]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+num_z]+=xp[n]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+num_z+N]+=-xp[n+N]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y-num_z]+=xp[n]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y-num_z+N]+=-xp[n+N]/(4*x_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+2*num_z]+=-xp[n+N]/(4*y_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-2*num_z]+=-xp[n+N]/(4*y_inc*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+2*num_z*num_y+N]+=xp[n]/(4*x_inc*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-2*num_z*num_y+N]+=xp[n]/(4*x_inc*x_inc)*2*DzetaDt[i][j][k]*C3;

			    //vertical advection

			    grad[n+2*N]+=((xp[n+num_z*num_y+1+N]-xp[n-num_z*num_y+1+N]-(xp[n+num_z*num_y-1+N]-xp[n-num_z*num_y-1+N]))/
				    (4*x_inc*z_inc)-((xp[n+num_z+1]-xp[n-num_z+1])-(xp[n+num_z-1]-xp[n-num_z-1]))/(4*y_inc*z_inc))
				*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+1+N]+=xp[n+2*N]/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+1+N]+=-xp[n+2*N]/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y-1+N]+=-xp[n+2*N]/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y-1+N]+=xp[n+2*N]/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z+1]+=-xp[n+2*N]/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z+1]+=xp[n+2*N]/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z-1]+=xp[n+2*N]/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z-1]+=-xp[n+2*N]/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;

			    //tilting terms

			    grad[n+num_z+2*N]+=-(xp[n+1]-xp[n-1])/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z+2*N]+=(xp[n+1]-xp[n-1])/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+1+N]+=(xp[n+num_z*num_y+2*N]-xp[n-num_z*num_y+2*N])/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-1+N]+=-(xp[n+num_z*num_y+2*N]-xp[n-num_z*num_y+2*N])/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y+2*N]+=(xp[n+1+N]-xp[n-1+N])/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+2*N]+=-(xp[n+1+N]-xp[n-1+N])/(4*x_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+1]+=-(xp[n+num_z+2*N]-xp[n-num_z+2*N])/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-1]+=(xp[n+num_z+2*N]-xp[n-num_z+2*N])/(4*y_inc*z_inc)*2*DzetaDt[i][j][k]*C3;

			    //stretching term (anelastic)
			    grad[n+num_z*num_y+N]+=((xp[n+num_z*num_y]-xp[n-num_z*num_y])/(2*x_inc)+
				    (xp[n+num_z+N]-xp[n-num_z+N])/(2*y_inc))/(2*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z*num_y+N]+=-((xp[n+num_z*num_y]-xp[n-num_z*num_y])/(2*x_inc)+
				    (xp[n+num_z+N]-xp[n-num_z+N])/(2*y_inc))/(2*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z]+=-((xp[n+num_z*num_y]-xp[n-num_z*num_y])/(2*x_inc)+
				    (xp[n+num_z+N]-xp[n-num_z+N])/(2*y_inc))/(2*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n-num_z]+=((xp[n+num_z*num_y]-xp[n-num_z*num_y])/(2*x_inc)+
				    (xp[n+num_z+N]-xp[n-num_z+N])/(2*y_inc))/(2*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y]+=((xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc)
				    -(xp[n+num_z]-xp[n-num_z])/(2*y_inc))/(2*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z*num_y]+=-((xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc)
				    -(xp[n+num_z]-xp[n-num_z])/(2*y_inc))/(2*x_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z+N]+=((xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc)
				    -(xp[n+num_z]-xp[n-num_z])/(2*y_inc))/(2*y_inc)*2*DzetaDt[i][j][k]*C3;
			    grad[n+num_z+N]+=-((xp[n+num_z*num_y+N]-xp[n-num_z*num_y+N])/(2*x_inc)
				    -(xp[n+num_z]-xp[n-num_z])/(2*y_inc))/(2*y_inc)*2*DzetaDt[i][j][k]*C3;
			}
		    }
		}
	    }
	}
    }

    // *** Adjoint for CalcSmooth

    if (C4b>0) {
	n=0;
	for (i=0; i < num_x; i++) {
	    for (j=0; j < num_y; j++) {
		for (k=0; k < num_z; k++) {
		    n++;
		    if (i>0&&i<=(num_dx-0)&&j>0&&j<=(num_dy-0)&&k>0&&k<=(num_dz-0)) {
			z=z_min+k*z_inc;
			if (z<=cutoff)
			    sss=ss;
			else
			    sss=1.0;

			if (z >= z_min) {
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C4=C4d;
			    else
				C4=C4a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C5=C5d;
			    else
				C5=C5a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C6=C6d;
			    else
				C6=C6a;
			    if (i<0||i>num_dx||j<0||j>num_dy||k<0||k>num_dz)
				C7=C7d;
			    else
				C7=C7a;
			    if (laplace==1 && i>0&&i<num_dx&&j>0&&j<num_dy&&k>0&&k<num_dz) {

				//Laplacian terms

				grad[n]+=(8*xp[n]-4*(xp[n+num_z*num_y]+xp[n-num_z*num_y]))/pow(x_inc,4)*C4*ss2;
				grad[n-num_z*num_y]+=(2*xp[n-num_z*num_y]-4*xp[n]+2*xp[n+num_z*num_y])/pow(x_inc,4)*C4*ss2; // dudx
				grad[n+num_z*num_y]+=(2*xp[n-num_z*num_y]-4*xp[n]+2*xp[n+num_z*num_y])/pow(x_inc,4)*C4*ss2;

				grad[n+N]+=(8*xp[n+N]-4*(xp[n+N+num_z*num_y]+xp[n+N-num_z*num_y]))/pow(x_inc,4)*C4*ss2;
				grad[n+N-num_z*num_y]+=(2*xp[n+N-num_z*num_y]-4*xp[n+N]+2*xp[n+N+num_z*num_y])/pow(x_inc,4)*C4*ss2; // dvdx
				grad[n+N+num_z*num_y]+=(2*xp[n+N-num_z*num_y]-4*xp[n+N]+2*xp[n+N+num_z*num_y])/pow(x_inc,4)*C4*ss2;

				grad[n+2*N]+=(8*xp[n+2*N]-4*(xp[n+2*N+num_z*num_y]+xp[n+2*N-num_z*num_y]))/pow(x_inc,4)*C7*sss;
				grad[n+2*N-num_z*num_y]+=(2*xp[n+2*N-num_z*num_y]-4*xp[n+2*N]+2*xp[n+2*N+num_z*num_y])/pow(x_inc,4)*C7*sss; // dwdx
				grad[n+2*N+num_z*num_y]+=(2*xp[n+2*N-num_z*num_y]-4*xp[n+2*N]+2*xp[n+2*N+num_z*num_y])/pow(x_inc,4)*C7*sss;

				grad[n]+=(8*xp[n]-4*(xp[n+num_z]+xp[n-num_z]))/pow(y_inc,4)*C4*ss2;
				grad[n-num_z]+=(2*xp[n-num_z]-4*xp[n]+2*xp[n+num_z])/pow(y_inc,4)*C4*ss2; // dudy
				grad[n+num_z]+=(2*xp[n-num_z]-4*xp[n]+2*xp[n+num_z])/pow(y_inc,4)*C4*ss2;

				grad[n+N]+=(8*xp[n+N]-4*(xp[n+N+num_z]+xp[n+N-num_z]))/pow(y_inc,4)*C4*ss2;
				grad[n+N-num_z]+=(2*xp[n+N-num_z]-4*xp[n+N]+2*xp[n+N+num_z])/pow(y_inc,4)*C4*ss2; // dvdy
				grad[n+N+num_z]+=(2*xp[n+N-num_z]-4*xp[n+N]+2*xp[n+N+num_z])/pow(y_inc,4)*C4*ss2;

				grad[n+2*N]+=(8*xp[n+2*N]-4*(xp[n+2*N+num_z]+xp[n+2*N-num_z]))/pow(y_inc,4)*C7*sss;
				grad[n+2*N-num_z]+=(2*xp[n+2*N-num_z]-4*xp[n+2*N]+2*xp[n+2*N+num_z])/pow(y_inc,4)*C7*sss; // dwdy
				grad[n+2*N+num_z]+=(2*xp[n+2*N-num_z]-4*xp[n+2*N]+2*xp[n+2*N+num_z])/pow(y_inc,4)*C7*sss;

				grad[n]+=(8*xp[n]-4*(xp[n+1]+xp[n-1]))/pow(z_inc,4)*C5*sss;
				grad[n-1]+=(2*xp[n-1]-4*xp[n]+2*xp[n+1])/pow(z_inc,4)*C5*sss; // dudz
				grad[n+1]+=(2*xp[n-1]-4*xp[n]+2*xp[n+1])/pow(z_inc,4)*C5*sss;

				grad[n+N]+=(8*xp[n+N]-4*(xp[n+N+1]+xp[n+N-1]))/pow(z_inc,4)*C5*sss;
				grad[n+N-1]+=(2*xp[n+N-1]-4*xp[n+N]+2*xp[n+N+1])/pow(z_inc,4)*C5*sss; // dvdz
				grad[n+N+1]+=(2*xp[n+N-1]-4*xp[n+N]+2*xp[n+N+1])/pow(z_inc,4)*C5*sss;

				grad[n+2*N]+=(8*xp[n+2*N]-4*(xp[n+2*N+1]+xp[n+2*N-1]))/pow(z_inc,4)*C6*sss;
				grad[n+2*N-1]+=(2*xp[n+2*N-1]-4*xp[n+2*N]+2*xp[n+2*N+1])/pow(z_inc,4)*C6*sss; // dwdz
				grad[n+2*N+1]+=(2*xp[n+2*N-1]-4*xp[n+2*N]+2*xp[n+2*N+1])/pow(z_inc,4)*C6*sss;

			    } // laplace=1

			    else if (laplace==0 && i>0 && j>0 && k>0) {

				//Gradient terms

				if (z<=cutoff) {
				    grad[n+2*N]+=2*(xp[n+2*N]-xp[n-1+2*N])/(z_inc*z_inc)*C6*ss;
				    grad[n-1+2*N]+=-2*(xp[n+2*N]-xp[n-1+2*N])/(z_inc*z_inc)*C6*ss;

				    grad[n]+=2*(xp[n]-xp[n-1])/(z_inc*z_inc)*C5*ss;
				    grad[n-1]+=-2*(xp[n]-xp[n-1])/(z_inc*z_inc)*C5*ss;
				    grad[n+N]+=2*(xp[n+N]-xp[n-1+N])/(z_inc*z_inc)*C5*ss;
				    grad[n-1+N]+=-2*(xp[n+N]-xp[n-1+N])/(z_inc*z_inc)*C5*ss;
				}

				else {
				    grad[n+2*N]+=2*(xp[n+2*N]-xp[n-1+2*N])/(z_inc*z_inc)*C6*ss;
				    grad[n-1+2*N]+=-2*(xp[n+2*N]-xp[n-1+2*N])/(z_inc*z_inc)*C6*ss;

				    grad[n]+=2*(xp[n]-xp[n-1])/(z_inc*z_inc)*C5*ss;
				    grad[n-1]+=-2*(xp[n]-xp[n-1])/(z_inc*z_inc)*C5*ss;
				    grad[n+N]+=2*(xp[n+N]-xp[n-1+N])/(z_inc*z_inc)*C5*ss;
				    grad[n-1+N]+=-2*(xp[n+N]-xp[n-1+N])/(z_inc*z_inc)*C5*ss;
				}

				grad[n+2*N]+=2*(xp[n+2*N]-xp[n-num_z*num_y+2*N])/(x_inc*x_inc)*C7*ss;
				grad[n-num_z*num_y+2*N]+=-2*(xp[n+2*N]-xp[n-num_z*num_y+2*N])/(x_inc*x_inc)*C7*ss;
				grad[n+2*N]+=2*(xp[n+2*N]-xp[n-num_z+2*N])/(y_inc*y_inc)*C7*ss;
				grad[n-num_z+2*N]+=-2*(xp[n+2*N]-xp[n-num_z+2*N])/(y_inc*y_inc)*C7*ss;

				grad[n]+=2*(xp[n]-xp[n-num_z*num_y])/(x_inc*x_inc)*C4*ss;
				grad[n-num_z*num_y]+=-2*(xp[n]-xp[n-num_z*num_y])/(x_inc*x_inc)*C4*ss;
				grad[n]+=2*(xp[n]-xp[n-num_z])/(y_inc*y_inc)*C4*ss;
				grad[n-num_z]+=-2*(xp[n]-xp[n-num_z])/(y_inc*y_inc)*C4*ss;

				grad[n+N]+=2*(xp[n+N]-xp[n-num_z*num_y+N])/(x_inc*x_inc)*C4*ss;
				grad[n-num_z*num_y+N]+=-2*(xp[n+N]-xp[n-num_z*num_y+N])/(x_inc*x_inc)*C4*ss;
				grad[n+N]+=2*(xp[n+N]-xp[n-num_z+N])/(y_inc*y_inc)*C4*ss;
				grad[n-num_z+N]+=-2*(xp[n+N]-xp[n-num_z+N])/(y_inc*y_inc)*C4*ss;

			    }
			}
		    }
		}
	    }
	}
    }

    // *** Adjoint for CalcBG

    if (C8b>0) {
	n=0;
	for (i=0; i < num_x; i++) {
	    for (j=0; j < num_y; j++) {
		for (k=0; k < num_z; k++) {
		    n++;
		    if ( i >= 2 && i < num_x-2
			    && j >= 2 && j < num_y-2
			    && k >= 0 && k < num_z ) {
			if ( coverage_bg[i][j][k] == 0 ) {
			    if (maxdbz[i][j][k]>10  && k<=weak_height){
				grad[n]+=2*C8a*(xp[n]-ub[n])/100.0;
				grad[n+N]+=2*C8a*(xp[n+N]-vb[n])/100.0;}
			    else {
				grad[n]+=2*C8a*(xp[n]-ub[n]);
				grad[n+N]+=2*C8a*(xp[n+N]-vb[n]);}
			}
		    }
		    if ( (i <= 1 || i >= num_x-1
				|| j <= 1 || j >= num_y-1)
			    && (k >= 0 && k < num_z) ) {
			if ( coverage_bg[i][j][k] == 0) {
			    if (maxdbz[i][j][k]>10  && k<=weak_height){
				grad[n]+=2*C8a*(xp[n]-ub[n])/100.0;
				grad[n+N]+=2*C8a*(xp[n+N]-vb[n])/100.0;}
			    else {
				grad[n]+=2*C8a*(xp[n]-ub[n]);
				grad[n+N]+=2*C8a*(xp[n+N]-vb[n]);}
			}

		    }
		}
	    }
	}
    }

    // impermeability condition
    n=0;
    for (int i=0; i < num_x; i++) {
	for (int j=0; j < num_y; j++) {
	    for (int k=0; k < num_z; k++) {
		n++;
		if (k==0)
		    grad[n+2*N]=0;
		if ((k==(num_z-1)) && upper_bc==1)
		    grad[n+2*N]=0;
	    }
	}
    }
}

void nrerror(const char *error_text)
    /* Numerical Recipes standard error handler */
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

void free_vector(double *v, long nl, long nh)
    /* free a double vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

double *vector(long nl, long nh)
{
    double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
}

// From Numerical Recipes in C++
double f1dim(double x)
{
    int j;
    double f,*xt;
    xt=vector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    f=(*nrfunc)(xt);
    free_vector(xt,1,ncom);
    return f;
}

// From Numerical Recipes in C++
double df1dim(double x)
{
    int j;
    double df1=0.0;
    double *xt,*df;
    xt=vector(1,ncom);
    df=vector(1,ncom);
    for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
    (*nrdfun)(xt,df);
    for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
    free_vector(df,1,ncom);
    free_vector(xt,1,ncom);
    return df1;
}

// From Numerical Recipes in C++
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
	double *fc, double (*func)(double))
{
    double ulim,u,r,q,fu,dum,largest;

    *fa=(*func)(*ax);
    *fb=(*func)(*bx);
    if (*fb > *fa) {
	SHFT(dum,*ax,*bx,dum)
	    SHFT(dum,*fb,*fa,dum)
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(*func)(*cx);
    while (*fb > *fc) {
	r=(*bx-*ax)*(*fb-*fc);
	q=(*bx-*cx)*(*fb-*fa);
	if (fabs(q-r)>TINY)
	    largest=fabs(q-r); // EDITTED
	else largest=TINY;
	u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(largest,q-r));
	ulim=(*bx)+GLIMIT*(*cx-*bx);
	if ((*bx-u)*(u-*cx) > 0.0) {
	    fu=(*func)(u);
	    if (fu < *fc) {
		*ax=(*bx);
		*bx=u;
		*fa=(*fb);
		*fb=fu;
		return;
	    } else if (fu > *fb) {
		*cx=u;
		*fc=fu;
		return;
	    }
	    u=(*cx)+GOLD*(*cx-*bx);
	    fu=(*func)(u);
	} else if ((*cx-u)*(u-ulim) > 0.0) {
	    fu=(*func)(u);
	    if (fu < *fc) {
		SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
		    SHFT(*fb,*fc,fu,(*func)(u))
	    }
	} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
	    u=ulim;
	    fu=(*func)(u);
	} else {
	    u=(*cx)+GOLD*(*cx-*bx);
	    fu=(*func)(u);
	}
	SHFT(*ax,*bx,*cx,u)
	    SHFT(*fa,*fb,*fc,fu)
    }
}

// From Numerical Recipes in C++
double dbrent(double ax, double bx, double cx, double (*f)(double),
	double (*df)(double), double tol, double *xmin)
{
    int iterb,ok1,ok2;
    double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
    double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

    a=(ax < cx ? ax : cx);
    b=(ax > cx ? ax : cx);
    x=w=v=bx;
    fw=fv=fx=(*f)(x);
    dw=dv=dx=(*df)(x); // FAILS
    for (iterb=1;iterb<=itmax_dbrent;iterb++) {
	xm=0.5*(a+b);
	tol1=tol*fabs(x)+ZEPS;
	tol2=2.0*tol1;
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
	    *xmin=x;
	    return fx;
	}
	if (fabs(e) > tol1) {
	    d1 = 2.0 * (b - a);
	    d = d2 = d1;
	    if (dw != dx) d1=(w-x)*dx/(dx-dw);
	    if (dv != dx) d2=(v-x)*dx/(dx-dv);
	    u1=x+d1;
	    u2=x+d2;
	    ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	    ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
	    olde=e;
	    e=d;
	    if (ok1 || ok2) {
		d = (ok1 && ok2) ? (fabs(d1) < fabs(d2) ? d1 : d2)
		    : (ok1) ? d1 : d2;
		if (fabs(d) <= fabs(0.5*olde)) {
		    u=x+d;
		    if (u-a < tol2 || b-u < tol2)
			d=SIGN(tol1,xm-x);
		} else {
		    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
	    } else {
		d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	    }
	} else {
	    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
	if (fabs(d) >= tol1) {
	    u=x+d;
	    fu=(*f)(u);
	} else {
	    u=x+SIGN(tol1,d);
	    fu=(*f)(u);
	    if (fu > fx) {
		*xmin=x;
		return fx;
	    }
	}
	du=(*df)(u);
	if (fu <= fx) {
	    if (u >= x) a=x; else b=x;
	    MOV3(v,fv,dv, w,fw,dw)
		MOV3(w,fw,dw, x,fx,dx)
		MOV3(x,fx,dx, u,fu,du)
	} else {
	    if (u < x) a=u; else b=u;
	    if (fu <= fw || w == x) {
		MOV3(v,fv,dv, w,fw,dw)
		    MOV3(w,fw,dw, u,fu,du)
	    } else if (fu < fv || v == x || v == w) {
		MOV3(v,fv,dv, u,fu,du)
	    }
	}
    }
    nrerror("Too many iterations in routine dbrent");
    return 0.0;
}

// From Numerical Recipes in C++
void dlinmin(double p[], double xi[], int n, double (*func)(double []),
	void (*dfunc)(double [], double []), double *mytemp_p)
{
    int j;
    double xx,xmin,fx,fb,fa,bx,ax;
    ncom=n;
    pcom=vector(1,n);
    xicom=vector(1,n);
    nrfunc=func;
    nrdfun=dfunc;
    for (j=1;j<=n;j++){
	pcom[j]=p[j];
	xicom[j]=xi[j];
    }

    ax=0.0;
    xx=1.0;
    mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);

    *mytemp_p = dbrent(ax, xx, bx, f1dim, df1dim, TOL, &xmin);

    for (j=1;j<=n;j++) {
	xi[j] *= xmin;
	p[j] += xi[j];
    }

    free_vector(xicom,1,n);
    free_vector(pcom,1,n);
}

// Adapted from Numerical Recipes in C++
void frprmn(double p[], int n, double ftol, int *iter, int itmax2,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
    double mytemp;
    int i, j, k, its;
    double gg,gam,dgg;
    double *g,*h,*xi;
    int nnn=0, stop=1;
    double max_diff;
    double mtoterr,norm_merr, totdenom;
    int a=0;
    double x, y, z,cov_count;
    double max_w;

    g=vector(1,n);
    h=vector(1,n);
    xi=vector(1,n);
    (*dfunc)(p,xi);
    for (j=1;j<=n;j++) {
	g[j] = -xi[j];
	xi[j]=h[j]=g[j];
    }
    for (i=0; i < num_x; i++) {
	for (j=0; j < num_y; j++) {
	    for (k=0; k < num_z; k++) {
		nnn++;
		p_old[nnn]=p[nnn]; // used in convergence criteria
	    }
	}
    }

    printf("      Its| MsCnsErr|  NrmzErr|AvgObsErr| MxObsErr|   ObsTot"
	    "|MsConsTot|  SmthTot|  VortTot|  MaxDiff|     MaxW");
    if ( vary_weights ) {
	printf("|      C2b|      C5b");
    }
    printf("\n");

    // minimization loop
    for (its = 1; its <= itmax2; its++) {
	*iter=its;
	dlinmin(p, xi, n, func, dfunc, &mytemp);
	output_error=1;
	CalcCost(p);
	output_error=1;
	nnn=0;
	stop=1;
	max_diff=0;
	cov_count=0;
	max_w=0.0;
	mtoterr = 0.0;
	totdenom = 0.0;
	for (k=0; k < (num_z); k++) {
	    for (i = 0; i < num_x; i++) {
		for (j = 0; j < num_y; j++) {
		    nnn++;
		    if ( coverage_bg[i][j][k] == 1 ) {
			mtoterr  += pow(div2[i][j][k],2);
			totdenom += dudx[i][j][k] * dudx[i][j][k]
			    + dvdy[i][j][k] * dvdy[i][j][k]
			    + dwdz[i][j][k] * dwdz[i][j][k];
			cov_count++;
			if ( p[nnn+2*N] > max_w ) {
			    max_w=p[nnn+2*N];
			}
		    }
		    if (its % 10 ==0) {
			if ( i > 1 && i < num_dx - 1
				&& j > 1 && j < num_dy - 1
				&& k >= 0 && k < num_z - 0 ) {

			    /*
			       Stop minimization prematurely if change in w
			       from 10 iterations ago is everywhere < 0.2 m/s
			     */

			    if (fabs(p[nnn + 2*N] - p_old[nnn + 2*N]) >= .2 ) {

			    }
			    if ( fabs(p[nnn+2*N]-p_old[nnn+2*N]) > max_diff ) {
				max_diff=fabs(p[nnn+2*N]-p_old[nnn+2*N]);
			    }

			}
			p_old[nnn+2*N] = p[nnn+2*N];
		    }
		}
	    }
	}
	if ((its % 10 == 0) && (its > 50) && (max_diff < 0.2)) {
	    stop = 0;
	}
	norm_merr=sqrt(mtoterr/totdenom);
	if (vary_weights==1){
	    if (its > 10 && meanobserror < 0.1) {
		if (norm_merr > mass_norm_thresh ){
		    C2a = C2a*1.2;
		    C2b = C2a;
		    C4a = C4a*1.2;
		    C3a = C3a*1.2;
		    C5a = C5a*1.2;
		    C3b = C3a;
		    C4b =C7a=C7b= C4a;
		    C5b =C6a=C6b= C5a;
		}
	    } // iteration if statement //
	    if (its > 100 && sqrt(maxobserror) > 5) {
		if (norm_merr < mass_norm_thresh ){
		    C2a = C2a*0.9;
		    C2b = C2a;
		    C4a = C4a*0.9;
		    C3a = C3a*0.9;
		    C5a = C5a*0.9;
		    C3b = C3a;
		    C4b =C7a=C7b= C4a;
		    C5b =C6a=C6b= C5a;
		}
	    } // iteration if statement //
	} // use adjustable weights if statement //

	if ( its % 10 == 0 ) {
	    printf("%9d|%9.4g|%9.4g|%9.4g|%9.4g|%9.4g|%9.4g|%9.4g|%9.4g|"
		    "%9.4g|%9.4g",
		    its, (C2b*mtoterr/cov_count), norm_merr,
		    meanobserror, sqrt(maxobserror),
		    obstotal/Jtotal, masstotal/Jtotal,
		    smoothtotal/Jtotal, vorttotal/Jtotal,
		    max_diff, max_w);
	    if ( vary_weights ) {
		printf("|%9.4g|%9.4g", C2b, C5b);
	    }
	    printf("\n");
	}
	if ( its > 9 && its % filt_freq == 0 && filt_type != NO_FILTER) {
	    switch (filt_type) {
		case NO_FILTER:
		    break;
		case LEISE:
		    Velfill(xp, num_x, num_y, num_z, 300, 10,1, 0.75);
		    VelLefilt(p, num_x, num_y, num_z, Lei_steps);
		    break;
		case LO_PASS:
		    Velfill(xp, num_x, num_y, num_z, 300, 10,1, 0.75);
		    Velfill(xp, num_x, num_y, num_z, 1, 1, 0, filt_alpha);
		    break;
	    }
	}
	norm_merr_old = norm_merr;
	meanobserror_old = meanobserror;

	if ( stop == 0 || its == itmax2 ) {
	    // Compute MAXDBZ field for netcdf output (untested code)
	    int NDIMS=3, NX=num_x, NY=num_dy+1, NZ=num_dz+1;
	    int ncid, x_dimid, y_dimid, z_dimid;
	    int uid, vid, wid, cvgid, cvgid2;
	    int xid, yid, zid, dbzid, retval;
	    int dimids[3], xdimid[1], ydimid[1], zdimid[1];
	    size_t index[3];
	    size_t ind;

	    printf("\nDONE \n");
	    if ((retval = nc_create(frprmn_fl_nm, NC_NETCDF4, &ncid)))
		ERR(retval);
	    if ((retval = nc_def_dim(ncid, "NX", NX, &x_dimid)))
		ERR(retval);
	    if ((retval = nc_def_dim(ncid, "NY", NY, &y_dimid)))
		ERR(retval);
	    if ((retval = nc_def_dim(ncid, "NZ", NZ, &z_dimid)))
		ERR(retval);
	    dimids[0] = x_dimid;
	    dimids[1] = y_dimid;
	    dimids[2] = z_dimid;
	    xdimid[0] = x_dimid;
	    ydimid[0] = y_dimid;
	    zdimid[0] = z_dimid;
	    if ((retval = nc_def_var(ncid, "X", NC_DOUBLE, 1,
			    xdimid, &xid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "Y", NC_DOUBLE, 1,
			    ydimid, &yid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "Z", NC_DOUBLE, 1,
			    zdimid, &zid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "U", NC_DOUBLE, NDIMS,
			    dimids, &uid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "V", NC_DOUBLE, NDIMS,
			    dimids, &vid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "W", NC_DOUBLE, NDIMS,
			    dimids, &wid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "CVG", NC_DOUBLE, NDIMS,
			    dimids, &cvgid)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "CVG2", NC_DOUBLE, NDIMS,
			    dimids, &cvgid2)))
		ERR(retval);
	    if ((retval = nc_def_var(ncid, "MAXDBZ", NC_DOUBLE, NDIMS,
			    dimids, &dbzid)))
		ERR(retval);
	    if ((retval = nc_enddef(ncid)))
		ERR(retval);

	    // Compute and output analyzed Vr (observation space)
	    Vr_error=1;
	    CalcCost(p);
	    Vr_error=0;
	    printf("stop=%d iterations=%d itmax2=%d \n", stop, its, itmax2);
	    for (k = 0; k < num_z; k++) {
		ind = k;
		z = z_min + k * z_inc;
		if ((retval = nc_put_var1_double(ncid, zid, &ind, &z)))
		    ERR(retval);
	    }
	    for (j=0; j < num_y; j++) {
		ind=j;
		y=y_min+j*y_inc;
		if ((retval = nc_put_var1_double(ncid, yid, &ind, &y)))
		    ERR(retval);
	    }
	    for (i=0; i < num_x; i++) {
		ind=i;
		x=x_min+i*x_inc;
		if ((retval = nc_put_var1_double(ncid, xid, &ind, &x)))
		    ERR(retval);
	    }
	    for (k=0; k<num_dz; k++) {
		for (j=0; j < num_y; j++) {
		    for (i=0; i < num_x; i++) {
			a=nn[i][j][k];
			index[0] = i;
			index[1] = j;
			index[2] = k;
			retval = nc_put_var1_double(ncid, uid, index, &p[a]);
			if ( retval != 0 ) {
			    ERR(retval);
			}
			retval = nc_put_var1_double(ncid, vid, index, &p[a+N]);
			if ( retval != 0 ) {
			    ERR(retval);
			}
			retval = nc_put_var1_double(ncid, wid, index,
				&p[a+2*N]);
			if ( retval != 0 ) {
			    ERR(retval);
			}
			retval = nc_put_var1_int(ncid, cvgid, index,
				&coverage_bg[i][j][k]);
			if ( retval != 0 ) {
			    ERR(retval);
			}
			retval = nc_put_var1_double(ncid, dbzid, index,
				&maxdbz[i][j][k]);
			if ( retval != 0 ) {
			    ERR(retval);
			}
		    }
		}
	    }
	    output_error=1;
	    if ((retval = nc_close(ncid)))
		ERR(retval);
	    printf("Convergence criteria met.\n");
	    printf(" ITS == %d\n", its);
	    return;
	}  //its if statement 5
	(*dfunc)(p,xi);

	// dgg += xi[j] * xi[j]; => Fletcher-Reeves
	// dgg += (xi[j] + g[j]) * xi[j]; => Polak-Ribiere
	dgg = gg = 0.0;
	for (j = 1; j <= n; j++) {
	    gg += g[j] * g[j];
	    dgg += (xi[j]+g[j])*xi[j];
	}
	if (gg == 0.0) {
	    printf("gradient is zero\n");
	    return;
	}
	gam=dgg/gg;
	for (j = 1; j <= n; j++) {
	    g[j] = -xi[j];
	    xi[j]=h[j]=g[j]+gam*h[j];
	}
    } // its minimization loop
    nrerror("Too many iterations in frprmn");
    printf("frprmn\n");
}

// Write out xp in netcdf
void writeout(double p[])
{
    double x,y,z;
    int NDIMS=3, NX=num_x, NY=num_dy+1, NZ=num_dz+1;
    int ncid, x_dimid, y_dimid, z_dimid;
    int uid, vid, wid, cvgid, cvgid2;
    int xid, yid, zid, dbzid, retval;
    int dimids[3], xdimid[1], ydimid[1], zdimid[1];
    size_t index[3];
    size_t ind;

    if ((retval = nc_create(writeout_fl_nm, NC_CLOBBER, &ncid)))
	ERR(retval);
    if ((retval = nc_def_dim(ncid, "NX", NX, &x_dimid)))
	ERR(retval);
    if ((retval = nc_def_dim(ncid, "NY", NY, &y_dimid)))
	ERR(retval);
    if ((retval = nc_def_dim(ncid, "NZ", NZ, &z_dimid)))
	ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    dimids[2] = z_dimid;
    xdimid[0] = x_dimid;
    ydimid[0] = y_dimid;
    zdimid[0] = z_dimid;

    if ((retval = nc_def_var(ncid, "X", NC_DOUBLE, 1,
		    xdimid, &xid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "Y", NC_DOUBLE, 1,
		    ydimid, &yid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "Z", NC_DOUBLE, 1,
		    zdimid, &zid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "U", NC_DOUBLE, NDIMS,
		    dimids, &uid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "V", NC_DOUBLE, NDIMS,
		    dimids, &vid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "W", NC_DOUBLE, NDIMS,
		    dimids, &wid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "CVG", NC_DOUBLE, NDIMS,
		    dimids, &cvgid)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "CVG2", NC_DOUBLE, NDIMS,
		    dimids, &cvgid2)))
	ERR(retval);
    if ((retval = nc_def_var(ncid, "MAXDBZ", NC_DOUBLE, NDIMS,
		    dimids, &dbzid)))
	ERR(retval);
    if ((retval = nc_enddef(ncid)))
	ERR(retval);


    int a=0;

    for (int k=0; k < num_z; k++) {
	ind=k;
	z =z_min + k * z_inc;
	if ((retval = nc_put_var1_double(ncid, zid, &ind, &z)))
	    ERR(retval);
    }// k loop
    for (int j=0; j < num_y; j++) {
	ind=j;
	y = y_min + j * y_inc;
	if ((retval = nc_put_var1_double(ncid, yid, &ind, &y)))
	    ERR(retval);
    } // j loop
    for (int i=0; i < num_x; i++) {
	ind=i;
	x = x_min + i * x_inc;
	if ((retval = nc_put_var1_double(ncid, xid, &ind, &x)))
	    ERR(retval);
    } // i loop

    for (int k=0; k<num_dz; k++) {
	for (int j=0; j < num_y; j++) {
	    for (int i=0; i < num_x; i++) {


		a=nn[i][j][k];

		index[0] = i;
		index[1] = j;
		index[2] = k;
		if ((retval = nc_put_var1_double(ncid, uid, index, &p[a])))
		    ERR(retval);
		if ((retval = nc_put_var1_double(ncid, vid, index, &p[a+N])))
		    ERR(retval);
		if ((retval = nc_put_var1_double(ncid, wid, index, &p[a+2*N])))
		    ERR(retval);
		if ((retval = nc_put_var1_int(ncid, cvgid, index, &coverage_bg[i][j][k])))
		    ERR(retval);
		if ((retval = nc_put_var1_double(ncid, dbzid, index, &maxdbz[i][j][k])))
		    ERR(retval);

	    }// i loop
	}// j loop
    }// k loop


    output_error=1;

    if ((retval = nc_close(ncid)))
	ERR(retval);

    return;

}

/*
   Adapted from code from Jidong Gao.

   Verify adjoint consistent with cost function - want to see fr decrease to
   1.0000 then zero
 */

void GradCheck(int n, double* xp)
{
    double *grad, *work;
    double rchek0=1e10, rchek,fx1,fx2,gxnn,ffff;
    int i,j;

    work = vector(1,n);
    grad = vector(1,n);
    rchek=rchek0;
    fx1 = CalcCost(xp);
    CalcGrad(xp,grad);
    rchek=rchek0;
    printf("GradCheck: rchek = %g. fx1 = %g\n", rchek, fx1);
    gxnn=0.0;
    for (j=1;j<=n;j++) {
	gxnn += grad[j]*grad[j];
    }
    printf("GradCheck: gxnn = %g\n", gxnn);
    for (j=1;j<=30;j++) {
	rchek=rchek*1e-1;
	for (i=1;i<=n;i++) {
	    work[i]=xp[i]+rchek*grad[i];
	}
	fx2 = CalcCost(work);
	ffff=(fx2-fx1)/(gxnn*rchek);
	printf("GradCheck: j = %d. fx2 = %.4g. ffff = %.4g\n", j, fx2, ffff);
    }
    printf("\n");
    printf("GradCheck: If above is all %f, trouble is coming.\n\n", NAN);
    free_vector(grad, 1, n);
    free_vector(work, 1, n);
}

void Filt2d(double **xy, double **fxy, int **miss,int nx, int ny, float alpha)
{
    int j, i;

    // copy the input field to the work array
    for (i=0; i<=nx; i++){
	for (j=0; j<=ny; j++) {
	    fxy[i][j] = xy[i][j];
	}
    }

    // filter the interior points
    for (i=2; i < nx; i++){
	for (j=2; j< ny; j++) {
	    if ( miss[i][j] == 0 || miss[i-1][j] == 0
		    || miss[i][j+1] == 0 || miss[i][j-1] == 0
		    || miss[i+1][j] == 0 || miss[i-1][j+1] == 0
		    || miss[i+1][j+1] == 0 || miss[i-1][j-1] == 0
		    || miss[i+1][j-1] == 0 ) {
		fxy[i][j] = xy[i][j] + 0.5 * alpha * (1.-alpha)
		    * (xy[i][j-1] + xy[i-1][j] + xy[i][j+1] + xy[i+1][j]
			    - 4.0 * xy[i][j])
		    + 0.25 * alpha * alpha
		    * (xy[i+1][j-1] + xy[i-1][j-1] + xy[i-1][j+1] + xy[i+1][j+1]
			    - 4.0 * xy[i][j]);
	    }
	}
    }

    // filter the x-boundaries
    for (i=1; i < nx; i++){
	fxy[i][0] = 0.5 * alpha * (xy[i+1][0] - 2.*xy[i][0] + xy[i-1][0])
	    + xy[i][0];
	fxy[i][ny] = 0.5 * alpha * (xy[i+1][ny] - 2.*xy[i][ny] + xy[i-1][ny])
	    + xy[i][ny];
    }

    // filter the y-boundaries
    for (j=1; j < ny; j++) {
	fxy[0][j] = 0.5 * alpha * (xy[0][j+1] - 2.*xy[0][j] + xy[0][j-1])
	    + xy[0][j];
	fxy[nx][j] = 0.5 * alpha * (xy[nx][j+1] - 2.*xy[nx][j] + xy[nx][j-1])
	    + xy[nx][j];
    }

    // set the corner values
    fxy[0][0]   = xy[0][0];
    fxy[0][ny]  = xy[0][ny];
    fxy[nx][0]  = xy[nx][0];
    fxy[nx][ny] = xy[nx][ny];

    for (i=0; i<=nx; i++){
	for (j=0; j<=ny; j++) {
	    xy[i][j] = fxy[i][j];
	}
    }
}

/*
   Horizontally filter missing velocities from 3DVAR wind analysis
 */
void Velfill(double *xp, int num_x, int num_y, int num_z, int uvit, int wit, int cov,float alpha)
{
    double **filtr;
    double **u2d;			/* u component */
    double **v2d;			/* v component */
    double **w2d;			/* w component */
    int **miss2d;
    //float alpha = 0.75;			/* lowpass smoothing parameter
    //				   (default = 0.5) */
    int a=0;				/* convert from 1D to 3D indices */
    int ifilt;
    int k, j, i;

    filtr = dalloc2(num_x, num_y, "filtr");
    u2d = dalloc2(num_x, num_y, "u2d");
    v2d = dalloc2(num_x, num_y, "v2d");
    w2d = dalloc2(num_x, num_y, "w2d");
    miss2d = ialloc2(num_x, num_y, "miss2d");

    for (k=0; k < num_z; k++) {

	/* Copy the flattened 1D radar fields into 2-D work arrays */
	for (i=0; i < num_x; i++){
	    for (j=0; j < num_y; j++) {
		a=nn[i][j][k];
		u2d[i][j] = xp[a];
		v2d[i][j] = xp[a+N];
		w2d[i][j] = xp[a+2*N];
		if (cov == 1){
		    miss2d[i][j] =  coverage_fil[i][j][k];}
		else miss2d[i][j] =  0;
	    }
	}

	// Filter missing data in 2-D u, v, and w arrays at level k
	for (ifilt = 1; ifilt <= uvit;ifilt++){
	    Filt2d(u2d,filtr,miss2d,num_dx,num_dy,alpha);
	    Filt2d(v2d,filtr,miss2d,num_dx,num_dy,alpha);
	}
	for (ifilt = 1; ifilt <= wit;ifilt++){
	    Filt2d(w2d,filtr,miss2d,num_dx,num_dy,alpha);
	}

	/* I (?) don't recommend filtering reflectivity where miss(i,j) = 1 */

	/* Copy 2-D work arrays back to 3-D field arrays */
	for (i=0; i < num_x; i++){
	    for (j=0; j < num_y; j++) {
		a=nn[i][j][k];
		xp[a] = u2d[i][j];
		xp[a+N] = v2d[i][j];
		xp[a+2*N] = w2d[i][j];

	    }
	}
    }

    dfree2(filtr);
    dfree2(u2d);
    dfree2(v2d);
    dfree2(w2d);
    ifree2(miss2d);
}
void VelLefilt(double *xp, int num_x, int num_y, int num_z, int steps)
{
    float u2dd[num_x][num_y];
    float v2dd[num_x][num_y];
    float w2dd[num_x][num_y];
    //float alpha = 0.75;			/* lowpass smoothing parameter
    //				   (default = 0.5) */
    int a=0;				/* convert from 1D to 3D indices */
    int k, j, i;
    int nz = 1;
    for (k=0; k < num_z; k++) {
	/* Copy the flattened 1D radar fields into 2-D work arrays */
	for (j=0; j < num_y; j++){
	    for (i=0; i < num_x; i++) {
		a=nn[i][j][k];
		u2dd[i][j] = xp[a];
		v2dd[i][j] = xp[a+N];
		w2dd[i][j] = xp[a+2*N];

	    } // j loop
	}  // i loop
	// Filter missing data in 2-D u, v, and w arrays at level k
	t5fltr_(u2dd[0], &num_y, &num_x, &nz,&steps);
	t5fltr_(v2dd[0], &num_y, &num_x, &nz,&steps);
	t5fltr_(w2dd[0], &num_y, &num_x, &nz,&steps);
	/* Copy 2-D work arrays back to 3-D field arrays */
	for (j=0; j < num_y; j++){
	    for (i=0; i < num_x; i++) {
		a=nn[i][j][k];
		xp[a]     = u2dd[i][j];
		xp[a+N]   = v2dd[i][j];
		xp[a+2*N] = w2dd[i][j];
	    }// j loop
	}// i loop

    }// k loop
    //dfree(u2d);
    //dfree(v2d);
    //dfree(w2d);
}
/*
   Return the distance between point (x1, y1) and point (x2, y2) on a plane.
 */

static double dist2(const double x1, const double y1, const double x2,
	const double y2)
{
    return hypot(x1 - x2, y1 - y2);
}

/*
   Allocate n int's. Initialize with 0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static int *ialloc(size_t n, const char *var_nm)
{
    int *arr;

    if ( !(arr = (int *)calloc(n, sizeof(int))) ) {
	fprintf(stderr, "Could not allocate %lu integer values for %s\n",
		(unsigned long)n, var_nm);
	exit(EXIT_FAILURE);
    }
    return arr;
}

/*
   Allocate n floats. Initialize with 0.0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static float *falloc(size_t n, const char *var_nm)
{
    float *arr;

    if ( !(arr = (float *)calloc(n, sizeof(float))) ) {
	fprintf(stderr, "Could not allocate %lu float values for %s\n",
		(unsigned long)n, var_nm);
	exit(EXIT_FAILURE);
    }
    return arr;
}

/*
   Allocate n doubles. Initialize with 0.0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static double *dalloc(size_t n, const char *var_nm)
{
    double *arr;

    if ( !(arr = (double *)calloc(n, sizeof(double))) ) {
	fprintf(stderr, "Could not allocate %lu double values for %s\n",
		(unsigned long)n, var_nm);
	exit(EXIT_FAILURE);
    }
    return arr;
}

/*
   Allocate allocate a two dimensional array of int's dimensioned
   [jmax][imax]. Initialize with INT_MIN. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static int **ialloc2(long jmax, long imax, const char *var_nm)
{
    int **dat = NULL;
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

    dat = (int **)calloc(jj + 2, sizeof(int **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (int *)calloc(jj * ii + 1, sizeof(int));
    if ( !dat[0] ) {
	fprintf(stderr, "Could not allocate 1st dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    for (j = 1; j <= jmax; j++) {
	dat[j] = dat[j - 1] + imax;
    }
    for (i = 0; i < jmax * imax; i++) {
	dat[0][i] = INT_MIN;
    }
    return dat;
}

/* Free memory allocated by dalloc2 */
static void ifree2(int **d)
{
    if ( d ) {
	free(d[0]);
    }
    free(d);
}

/*
   Allocate allocate a two dimensional array of double's dimensioned
   [jmax][imax]. Initialize with 0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static double **dalloc2(long jmax, long imax, const char *var_nm)
{
    double **dat = NULL;
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

    dat = (double **)calloc(jj + 2, sizeof(double **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (double *)calloc(jj * ii + 1, sizeof(double));
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

/* Free memory allocated by dalloc2 */
static void dfree2(double **d)
{
    if ( d ) {
	free(d[0]);
    }
    free(d);
}

/*
   Allocate allocate a four dimensional array of int's dimensioned
   [kmax][jmax][imax]. Initialize with 0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static int ***ialloc3(long kmax, long jmax, long imax, const char *var_nm)
{
    int ***dat = NULL;
    long k, j, i;
    size_t kk, jj, ii;

    /* Make sure casting to size_t does not overflow anything.  */
    if (kmax <= 0 || jmax <= 0 || imax <= 0) {
	fprintf(stderr, "Array dimensions must be positive for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    kk = (size_t)kmax;
    jj = (size_t)jmax;
    ii = (size_t)imax;
    if ((kk * jj) / kk != jj || (kk * jj * ii) / (kk * jj) != ii) {
	fprintf(stderr, "Dimensions [%ld][%ld][%ld]  too big "
		"for pointer arithmetic when allocating %s.\n",
		kmax, jmax, imax, var_nm);
	exit(EXIT_FAILURE);
    }

    dat = (int ***)calloc(kk + 2, sizeof(int **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (int **)calloc(kk * jj + 1, sizeof(int *));
    if ( !dat[0] ) {
	fprintf(stderr, "Could not allocate 1st dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0][0] = (int *)calloc(kk * jj * ii, sizeof(int));
    if ( !dat[0][0] ) {
	fprintf(stderr, "Could not allocate array of values for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    for (k = 1; k <= kmax; k++) {
	dat[k] = dat[k - 1] + jmax;
    }
    for (j = 1; j <= kmax * jmax; j++) {
	dat[0][j] = dat[0][j - 1] + imax;
    }
    for (i = 0; i < kmax * jmax * imax; i++) {
	dat[0][0][i] = 0;
    }
    return dat;
}

/*
   Allocate allocate a four dimensional array of float's dimensioned
   [kmax][jmax][imax]. Initialize with 0.0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static float ***falloc3(long kmax, long jmax, long imax, const char *var_nm)
{
    float ***dat = NULL;
    long k, j, i;
    size_t kk, jj, ii;

    /* Make sure casting to size_t does not overflow anything.  */
    if (kmax <= 0 || jmax <= 0 || imax <= 0) {
	fprintf(stderr, "Array dimensions must be positive for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    kk = (size_t)kmax;
    jj = (size_t)jmax;
    ii = (size_t)imax;
    if ((kk * jj) / kk != jj || (kk * jj * ii) / (kk * jj) != ii) {
	fprintf(stderr, "Dimensions [%ld][%ld][%ld]  too big "
		"for pointer arithmetic when allocating %s.\n",
		kmax, jmax, imax, var_nm);
	exit(EXIT_FAILURE);
    }

    dat = (float ***)calloc(kk + 2, sizeof(float **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (float **)calloc(kk * jj + 1, sizeof(float *));
    if ( !dat[0] ) {
	fprintf(stderr, "Could not allocate 1st dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0][0] = (float *)calloc(kk * jj * ii, sizeof(float));
    if ( !dat[0][0] ) {
	fprintf(stderr, "Could not allocate array of values for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    for (k = 1; k <= kmax; k++) {
	dat[k] = dat[k - 1] + jmax;
    }
    for (j = 1; j <= kmax * jmax; j++) {
	dat[0][j] = dat[0][j - 1] + imax;
    }
    for (i = 0; i < kmax * jmax * imax; i++) {
	dat[0][0][i] = 0.0;
    }
    return dat;
}

/*
   Allocate allocate a three dimensional array of double's dimensioned
   [kmax][jmax][imax]. Initialize with 0.0. Exit on failure.
   var_nm, for error message, gives name of variable receiving allocation.
 */

static double ***dalloc3(long kmax, long jmax, long imax, const char *var_nm)
{
    double ***dat = NULL;
    long k, j, i;
    size_t kk, jj, ii;

    /* Make sure casting to size_t does not overflow anything.  */
    if (kmax <= 0 || jmax <= 0 || imax <= 0) {
	fprintf(stderr, "Array dimensions must be positive for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    kk = (size_t)kmax;
    jj = (size_t)jmax;
    ii = (size_t)imax;
    if ((kk * jj) / kk != jj || (kk * jj * ii) / (kk * jj) != ii) {
	fprintf(stderr, "Dimensions [%ld][%ld][%ld]  too big "
		"for pointer arithmetic when allocating %s.\n",
		kmax, jmax, imax, var_nm);
	exit(EXIT_FAILURE);
    }

    dat = (double ***)calloc(kk + 2, sizeof(double **));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 2nd dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (double **)calloc(kk * jj + 1, sizeof(double *));
    if ( !dat[0] ) {
	fprintf(stderr, "Could not allocate 1st dimension for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    dat[0][0] = (double *)calloc(kk * jj * ii, sizeof(double));
    if ( !dat[0][0] ) {
	fprintf(stderr, "Could not allocate array of values for %s.\n",
		var_nm);
	exit(EXIT_FAILURE);
    }
    for (k = 1; k <= kmax; k++) {
	dat[k] = dat[k - 1] + jmax;
    }
    for (j = 1; j <= kmax * jmax; j++) {
	dat[0][j] = dat[0][j - 1] + imax;
    }
    for (i = 0; i < kmax * jmax * imax; i++) {
	dat[0][0][i] = 0.0;
    }
    return dat;
}

/*
   Allocate an array dimensioned [lmax][kmax][jmax][imax].
   Only allocate the first, second, and third dimensions. The fourth
   dimension is allocated as d. In other words, if return value is dat,
   dat[0][0][0] = d.  Array is named nm in error messages.
   Exit if anything goes wrong.
 */

static float ****calloc4f(long lmax, long kmax, long jmax, long imax, float *d,
	char *nm)
{
    float ****dat;
    long k, j, l;
    size_t ll, kk, jj, ii;

    /* Make sure casting to size_t does not overflow anything.  */
    if (lmax <= 0 || kmax <= 0 || jmax <= 0 || imax <= 0) {
	fprintf(stderr, "Array dimensions must be positive in %s.\n", nm);
	exit(EXIT_FAILURE);
    }
    ll = (size_t)lmax;
    kk = (size_t)kmax;
    jj = (size_t)jmax;
    ii = (size_t)imax;
    if ((ll * kk) / ll != kk || (ll * kk * jj) / (ll * kk) != jj
	    || (ll * kk * jj * ii) / (ll * kk * jj) != ii) {
	fprintf(stderr, "Dimensions [%ld][%ld][%ld][%ld]  too big "
		"for pointer arithmetic in %s.\n", lmax, kmax, jmax, imax, nm);
	exit(EXIT_FAILURE);
    }

    dat = (float ****)CALLOC(ll + 2, sizeof(float ***));
    if ( !dat ) {
	fprintf(stderr, "Could not allocate 3rd dimension of %s.\n", nm);
	exit(EXIT_FAILURE);
    }
    dat[0] = (float ***)CALLOC(ll * kk + 1, sizeof(float **));
    if ( !dat[0] ) {
	FREE(dat);
	fprintf(stderr, "Could not allocate 2nd dimension of %s.\n", nm);
	exit(EXIT_FAILURE);
    }
    dat[0][0] = (float **)CALLOC(ll * kk * jj + 1, sizeof(float *));
    if ( !dat[0][0] ) {
	FREE(dat[0]);
	FREE(dat);
	fprintf(stderr, "Could not allocate 1st dimension of %s.\n", nm);
	exit(EXIT_FAILURE);
    }
    dat[0][0][0] = d;
    for (l = 1; l <= lmax; l++) {
	dat[l] = dat[l - 1] + kmax;
    }
    for (k = 1; k <= lmax * kmax; k++) {
	dat[0][k] = dat[0][k - 1] + jmax;
    }
    for (j = 1; j <= lmax * kmax * jmax; j++) {
	dat[0][0][j] = dat[0][0][j - 1] + imax;
    }
    return dat;
}

/* Free memory allocated with alloc4f */
static void free4f(float ****dat)
{
    if (dat) {
	if (dat[0]) {
	    if (dat[0][0]) {
		FREE(dat[0][0]);
	    }
	    FREE(dat[0]);
	}
	FREE(dat);
    }
}
