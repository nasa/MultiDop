/*
   -	geog_lib.h --
   -		Declarations of structures and functions
   -		that store and manipulate geographic data.
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
   .	$Revision: 1.29 $ $Date: 2013/05/13 22:54:16 $
 */

#ifndef GEOG_LIB_H_
#define GEOG_LIB_H_

#define GEOG_VERSION "1.0"

#include <stdlib.h>

#ifndef M_PI
#define M_PI     3.141592653589793238462	/* pi */
#endif
#ifndef M_PI_2
#define M_PI_2   1.570796326794896619231	/* pi / 2 */
#endif
#ifndef M_PI_4
#define M_PI_4   0.78539816339744830961		/* pi / 4 */
#endif
#ifndef RAD_DEG
#define RAD_DEG   0.01745329251994329576	/* radians / degree */
#endif
#ifndef DEG_RAD
#define DEG_RAD   57.29577951308232087680	/* degrees / radian */
#endif

/* A geographic point */ 
struct GeogPt {
    double lon;			/* Longitude, radians */
    double lat;			/* Latitude, radians */
};

void GeogDMS(double, double *, double  *, double *, char *);
double GeogREarth(const double *);
double GeogLonR(const double, const double);
double GeogLonDiff(const double, const double);
double GeogLatN(const double);
double GeogDist(const double, const double, const double, const double);
double GeogAz(const double, const double, const double, const double);
void GeogStep(const double, const double, const double, const double,
	double *, double *);
double GeogBeamHt(double, double, double);
int GeogContainPt(const struct GeogPt, const struct GeogPt *, const size_t);

#endif
