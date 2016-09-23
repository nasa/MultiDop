/*
   -	geog_lib.c --
   -		This file defines functions that do geography
   -		calculations and manage geographic data.
   -		See geog_lib (3).
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
   .	$Revision: 1.39 $ $Date: 2013/05/13 22:54:16 $
   .
   .	References
   .	Smart, W. M., "Textbook on Spherical Astronomy",
   .	Sixth edition revised by R. M. Green.
   .	Cambridge University Press, Cambridge, 1977.
   .
   .	Sinnott, R. W., "Virtues of the Haversine",
   .	Sky and Telescope, vol. 68, no. 2, 1984, p. 159
   .	cited in: http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include "geog_lib.h"

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

/* Convert decimal degrees to degrees-minutes-seconds */
void GeogDMS(double ddeg, double *deg, double *min, double *sec, char *fmt)
{
    double m;			/* Minutes, including fraction */
    char buf[LDBL_MAX_10_EXP];	/* Character representations of minutes or
				   seconds, formatted with fmt */
    double v;			/* Minutes or seconds, read from buf */

    m = modf(ddeg, deg) * 60.0;
    *sec = modf(m, min) * 60.0;
    if ( fmt ) {
	snprintf(buf, LDBL_MAX_10_EXP, fmt, *sec);
	v = strtod(buf, NULL);
	if ( v == 60.0 ) {
	    *sec = 0.0;
	    *min += 1.0;
	}
	if ( v == -60.0 ) {
	    *sec = 0.0;
	    *min -= 1.0;
	}
	snprintf(buf, LDBL_MAX_10_EXP, fmt, *min);
	v = strtod(buf, NULL);
	if ( v == 60.0 ) {
	    *min = 0.0;
	    *deg += 1.0;
	}
	if ( v == -60.0 ) {
	    *min = 0.0;
	    *deg -= 1.0;
	}
    }
}

/* Get or set Earth radius */
double GeogREarth(const double *r)
{
    static double rearth = 6366707.019;		/* 1' = 1852 m */
    static int init = 0;
    char *s;

    if ( r ) {
	rearth = *r;
	init = 1;
    } else if ( !init && (s = getenv("GEOG_REARTH")) ) {
	rearth = strtod(s, NULL);
	init = 1;
    }
    return rearth;
}

/* Put l into [r - M_PI, r + M_PI) */
double GeogLonR(const double l, const double r)
{
    double l1 = fmod(l, 2.0 * M_PI);
    l1 = (l1 < r - M_PI) ? l1 + 2.0 * M_PI
	: (l1 >= r + M_PI) ? l1 - 2.0 * M_PI : l1;
    return (l1 == -0.0) ? 0.0 : l1;
}

/* Return l1 - l0 in [l0 - M_PI, l0 + M_PI) */
double GeogLonDiff(const double l1, const double l0)
{
    return GeogLonR(l1, l0) - l0;
}

/* Go l radians north of equator */
double GeogLatN(const double l)
{
    double l1 = fmod(l, 2.0 * M_PI);
    l1 += (l1 < 0.0) ? 2.0 * M_PI : 0.0;
    return (l1 > 1.5 * M_PI) ? l1 - 2.0 * M_PI
	: (l1 > M_PI_2 ) ? M_PI - l1 : l1;
}

/* Great circle distance in radians between two points */
double GeogDist(const double o1, const double a1, const double o2,
	const double a2)
{
    double sin_do_2, sin_da_2, a;

    sin_do_2 = sin(0.5 * (o2 - o1));
    sin_da_2 = sin(0.5 * (a2 - a1));
    a = sqrt(sin_da_2 * sin_da_2 + cos(a1) * cos(a2) * sin_do_2 * sin_do_2);
    return (a > 1.0 ? M_PI : 2.0 * asin(a));
}

/* GeogAz from (longitude, latitude): (o1, a1) to (o2, a2) */
double GeogAz(const double o1, const double a1, const double o2,
	const double a2)
{
    double sin_da, sin_sa, y, x;

    sin_da = sin(a1 - a2);
    sin_sa = sin(a2 + a1);
    y = cos(a2) * sin(o2 - o1);
    x = 0.5 * (sin_sa - sin_da - (sin_sa + sin_da) * cos(o2 - o1));
    return atan2(y, x);
}
/*
   Compute destination point longitude *o2, latitude *a2 at given separation s
   and direction d from point at longitude = o1, latitude a1.
 */
void GeogStep(const double o1, const double a1, const double d, const double s,
	double *o2, double *a2)
{
    double sin_s, sin_d, cos_d, dlon, a, x, y;

    sin_s = sin(s);
    sin_d = sin(d);
    cos_d = cos(d);
    a = 0.5 * (sin(a1 + s) * (1.0 + cos_d) + sin(a1 - s) * (1.0 - cos_d));
    *a2 = (a > 1.0) ? M_PI_2 : (a < -1.0) ? -M_PI_2 : asin(a);
    y = sin_s * sin_d;
    x = 0.5 * (cos(a1 + s) * (1 + cos_d) + cos(a1 - s) * (1 - cos_d));
    dlon = atan2(y, x);
    *o2 = GeogLonR(o1 + dlon, 0.0);
}

/*
   Height above ground after traveling distance d along a line tilt radians
   above horizontal. a0 is radius of Earth. a0 and d must have same units, which
   will be unit of answer.
 */

double GeogBeamHt(double d, double tilt, double a0)
{
    return sqrt(a0 * a0 + 2 * a0 * d * sin(tilt) + d * d) - a0;
}

/*
   This function returns true if polygon pts contains point at longitude lon,
   latitude lat.  pts is an array of lon0, lat0, lon1, lat1, ... values for
   n_pts points.
 */

int GeogContainPt(const struct GeogPt pt, const struct GeogPt *pts,
	const size_t n_pts)
{
    int mrdx;				/* Number of times a line crosses
					   meridian containing (lon, lat) */
    int lnx;				/* Number of times a line crosses line
					   from (lon, lat) to North pole */
    const struct GeogPt *p0, *p1;	/* Points from pts */
    double lon0, lat0;			/* p0 */
    double lon1, lat1;			/* p1 */
    double z;				/* Distance along Earth's axis, from
					   center of Earth */

    /*
       Loop through segments in pts, counting number of times pts
       crosses meridian, and number of crossings between (lon, lat) and
       North pole
     */

    for (mrdx = lnx = 0, p0 = pts + n_pts - 1, p1 = pts;
	    p1 < pts + n_pts; p0 = p1++) {

	/*
	   Determine if segment defined by p0--p1 straddles meridian
	   containing geoPt, or is on boundary.  Do not count segments on
	   boundary as more than one crossing
	 */

	lon0 = GeogLonR(p0->lon, pt.lon);
	lon1 = GeogLonR(p1->lon, pt.lon);
	if ( ( fabs(lon0 - lon1) < M_PI
		    && (   (lon0 < pt.lon && pt.lon <= lon1)
			|| (lon1 < pt.lon && pt.lon <= lon0))) ) {
	    double xlat;		/* Latitude of segment crossing */

	    mrdx++;
	    lat0 = p0->lat;
	    lat1 = p1->lat;
	    xlat = lat0 + (pt.lon - lon0) * (lat1 - lat0) / (lon1 - lon0);
	    if ( xlat > pt.lat ) {
		lnx = !lnx;
	    }

	}
    }

    if ( mrdx % 2 == 1 ) {
	/*
	   Odd number of meridian crossings => region contains a pole.
	   Assume pole is that of hemisphere containing pts mean.
	   If pts polygon contains North Pole, negate the result.
	 */

	for (p0 = pts, z = 0.0; p0 < pts + n_pts; p0++) {
	    z += sin(p0->lat);
	}
	if ( z > 0.0 ) {
	    lnx = !lnx;
	}

    }
    return lnx;
}

