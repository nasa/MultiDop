/*
   -	geog_proj.c --
   -		This file defines functions that convert back and forth between
   -		longitude latitude pairs and map points.  See geog_proj (3).
   - 
   .	Copyright (c) 2012 Gordon D. Carrie.  All rights reserved.
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
   .	$Revision: 1.6 $ $Date: 2013/05/10 22:32:05 $
 */

/*

   Ref.
   Snyder, John P.
   Map Projections used by the U.S. Geological Survey.
   (Geological Survey bulletin ; 1532)
   United States Government Printing Office, Washington:  1982.

 */

#include <math.h>
#include <stdio.h>
#include "geog_lib.h"
#include "geog_proj.h"

static struct GeogProj setRefPtProj(double, double);

int GeogProjLonLatToXY(double lon, double lat, double *x_p, double *y_p,
	struct GeogProj *projPtr)
{
    switch (projPtr->type) {
	case CylEqDist:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;

		*x_p = GeogLonDiff(lon, lon0) * cos_lat0 * r0;
		*y_p = (lat - lat0) * r0;
	    }
	    break;
	case CylEqArea:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.lon0;

		*x_p = r0 * GeogLonDiff(lon, lon0);
		*y_p = r0 * sin(lat);
	    }
	    break;
	case Mercator:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.lon0;
		double limit;

		limit = M_PI_2 * 8.0 / 9.0;	/* 80 degrees */
		if ( fabs(lat) > limit ) {
		    *x_p = *y_p = NAN;
		    return 0;
		}
		*x_p = r0 * GeogLonDiff(lon, lon0);
		*y_p = r0 * log(tan(M_PI_4 + 0.5 * lat));

	    }
	    break;
	case LambertConfConic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.LambertConfConic.lon0;
		double F = projPtr->params.LambertConfConic.F;
		double n = projPtr->params.LambertConfConic.n;
		double rho0 = projPtr->params.LambertConfConic.rho0;
		double rho, theta;

		rho = r0 * F / pow(tan(M_PI_4 + 0.5 * lat), n);
		theta = n * GeogLonDiff(lon, lon0);
		*x_p = rho * sin(theta);
		*y_p = rho0 - rho * cos(theta);
	    }
	    break;
	case LambertEqArea:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double k, dlon;
		double cos_lat, sin_lat, cos_dlon;

		cos_lat = cos(lat);
		sin_lat = sin(lat);
		lon = GeogLonR(lon, lon0);
		dlon = lon - lon0;
		cos_dlon = cos(dlon);
		if ( GeogDist(lon0, lat0, lon, lat) > M_PI_2 ) {
		    *x_p = *y_p = NAN;
		    return 0;
		}
		k = 1.0 + sin_lat0 * sin_lat + cos_lat0 * cos_lat * cos_dlon;
		k = sqrt(2.0 / k);
		*x_p = r0 * k * cos_lat * sin(dlon);
		*y_p = r0 * k
		    * (cos_lat0 * sin_lat - sin_lat0 * cos_lat * cos_dlon);
	    }
	    break;
	case Orthographic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double cos_lat, dlon;

		if ( GeogDist(lon0, lat0, lon, lat) > M_PI_2 ) {
		    *x_p = *y_p = NAN;
		    return 0;
		}
		cos_lat = cos(lat);
		dlon =  GeogLonDiff(lon, lon0);
		*x_p = r0 * cos_lat * sin(dlon);
		*y_p = r0 * (cos_lat0 * sin(lat)
			- sin_lat0 * cos_lat * cos(dlon));
	    }
	    break;
	case Stereographic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double dlon, cos_dlon, k, cos_lat, sin_lat;

		/*
		   Follow convention and treat as hemisphere projection.
		 */

		if ( GeogDist(lon0, lat0, lon, lat) > M_PI_2 ) {
		    *x_p = *y_p = NAN;
		    return 0;
		}
		cos_lat = cos(lat);
		sin_lat = sin(lat);
		dlon = GeogLonDiff(lon, lon0);
		cos_dlon = cos(dlon);
		k = 2.0 / (1.0 + sin_lat0 * sin_lat
			+ cos_lat0 * cos_lat * cos_dlon);
		*x_p = r0 * k * cos_lat * sin(dlon);
		*y_p = r0 * k
		    * (cos_lat0 * sin_lat - sin_lat0 * cos_lat * cos_dlon);
	    }
	    break;
    }
    if (projPtr->rotation != 0) {
	double x = *x_p, y = *y_p;
	double x_;

	x_ = x * projPtr->cosr + y * projPtr->sinr;
	*y_p = y * projPtr->cosr - x * projPtr->sinr;
	*x_p = x_;
    }
    return 1;
}

int GeogProjXYToLonLat(double x, double y, double *lon_p, double *lat_p,
	struct GeogProj *projPtr)
{
    if (projPtr->rotation != 0) {
	double x_;

	x_ =  x * projPtr->cosr - y * projPtr->sinr;
	y = x * projPtr->sinr + y * projPtr->cosr;
	x = x_;
    }
    switch (projPtr->type) {
	case CylEqDist:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;

		*lon_p = GeogLonR(lon0 + x / (cos_lat0 * r0), lon0);
		*lat_p = y / r0;
	    }
	    break;
	case CylEqArea:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.lon0;
		double r;

		r = GeogLatN(y / r0);
		*lat_p = asin(r);
		*lon_p = GeogLonR(lon0 + x / r0, lon0);
	    }
	    break;
	case Mercator:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.lon0;

		*lon_p = GeogLonR(lon0 + x / r0, lon0);
		*lat_p = M_PI_2 - 2.0 * atan(exp(-y / r0));
	    }
	    break;
	case LambertConfConic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.LambertConfConic.lon0;
		double F = projPtr->params.LambertConfConic.F;
		double n = projPtr->params.LambertConfConic.n;
		double rho0 = projPtr->params.LambertConfConic.rho0;
		double rho, theta;

		rho = hypot(x, rho0 - y);
		rho = copysign(rho, n);
		if ( n < 0.0 ) {
		    x = -x;
		    y = -y;
		    rho0 = -rho0;
		}
		theta = atan2(x, rho0 - y);
		*lon_p = GeogLonR(theta / n + lon0, lon0);
		if ( rho != 0.0 ) {
		    *lat_p = 2.0 * atan(pow(r0 * F / rho, 1.0 / n)) - M_PI_2;
		} else {
		    *lat_p = M_PI_2;
		    *lat_p = copysign(*lat_p, n);
		}
	    }
	    break;
	case LambertEqArea:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double rho, c, cos_c, sin_c, ord, lon, lat;

		rho = hypot(x, y);
		if (rho > 2.0 * r0) {
		    *lon_p = *lat_p = NAN;
		    return 0;
		}
		c = 2.0 * asin(rho / (2.0 * r0));
		cos_c = cos(c);
		sin_c = sin(c);
		if (rho == 0.0) {
		    lat = lat0;
		} else {
		    ord = cos_c * sin_lat0 + (y * sin_c * cos_lat0 / rho);
		    if (ord > 1.0)  {
			*lon_p = *lat_p = NAN;
			return 0;
		    }
		    lat = asin(ord);
		}
		*lat_p = lat;
		lon = lon0 + atan2(x * sin_c,
			(rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c));
		*lon_p = GeogLonR(lon, lon0);
	    }
	    break;
	case Orthographic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double rho, c, cos_c, sin_c, ord, lon, lat;

		rho = hypot(x, y);
		if ( rho / r0 > 1.0 ) {
		    *lon_p = *lat_p = NAN;
		    return 0;
		}
		c = asin(rho / r0);
		cos_c = cos(c);
		sin_c = sin(c);
		if ( rho == 0.0 ) {
		    lat = lat0;
		} else {
		    ord = cos_c * sin_lat0 + (y * sin_c * cos_lat0 / rho);
		    if ( ord > 1.0 ) {
			*lon_p = *lat_p = NAN;
			return 0;
		    }
		    lat = asin(ord);
		}
		*lat_p = lat;
		lon = lon0 + atan2(x * sin_c,
			(rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c));
		*lon_p = GeogLonR(lon, lon0);
	    }
	    break;
	case Stereographic:
	    {
		double r0 = GeogREarth(NULL);
		double lon0 = projPtr->params.RefPt.lon0;
		double lat0 = projPtr->params.RefPt.lat0;
		double cos_lat0 = projPtr->params.RefPt.cos_lat0;
		double sin_lat0 = projPtr->params.RefPt.sin_lat0;
		double rho, c, cos_c, sin_c, ord, lon, lat;

		rho = hypot(x, y);
		c = 2.0 * atan2(rho, 2.0 * r0);
		cos_c = cos(c);
		sin_c = sin(c);
		if (rho == 0.0) {
		    lat = lat0;
		} else {
		    ord = cos_c * sin_lat0 + (y * sin_c * cos_lat0 / rho);
		    if (ord > 1.0) {
			*lon_p = *lat_p = NAN;
			return 0;
		    }
		    lat = asin(ord);
		}
		*lat_p = lat;
		lon = lon0 + atan2(x * sin_c,
			rho * cos_lat0 * cos_c - y * sin_lat0 * sin_c);
		*lon_p = GeogLonR(lon, lon0);
	    }
	    break;
    }
    return 1;
}

void GeogProjSetRotation(struct GeogProj *projPtr, double angle)
{
    projPtr->rotation = angle;
    projPtr->cosr = cos(angle);
    projPtr->sinr = sin(angle);
}

int GeogProjSetFmStr(char *l, struct GeogProj *projPtr)
{
    double lon0, lat0, lat1, lat2;

    if ( sscanf(l, "CylEqDist %lf %lf", &lon0, &lat0) == 2 ) {
	return GeogProjSetCylEqDist(lon0 * RAD_DEG, lat0 * RAD_DEG, projPtr);
    } else if ( sscanf(l, "CylEqArea %lf", &lon0) == 1 ) {
	return GeogProjSetCylEqArea(lon0 * RAD_DEG, projPtr);
    } else if ( sscanf(l, "Mercator %lf", &lon0) == 1 ) {
	return GeogProjSetMercator(lon0 * RAD_DEG, projPtr);
    } else if ( sscanf(l, "LambertConfConic %lf %lf %lf %lf",
		&lon0, &lat0, &lat1, &lat2) == 4 ) {
	return GeogProjSetLambertConfConic(lon0 * RAD_DEG, lat0 * RAD_DEG,
		lat1 * RAD_DEG, lat2 * RAD_DEG, projPtr);
    } else if ( sscanf(l, "LambertEqArea %lf %lf", &lon0, &lat0) == 2 ) {
	return GeogProjSetLambertEqArea(lon0 * RAD_DEG, lat0 * RAD_DEG,
		projPtr);
    } else if ( sscanf(l, "Stereographic %lf %lf", &lon0, &lat0) == 2 ) {
	return GeogProjSetStereographic(lon0 * RAD_DEG, lat0 * RAD_DEG,
		projPtr);
    } else if ( sscanf(l, "Orthographic %lf %lf", &lon0, &lat0) == 2 ) {
	return GeogProjSetOrthographic(lon0 * RAD_DEG, lat0 * RAD_DEG, projPtr);
    }
    return 0;
}

int GeogProjSetCylEqDist(double lon0, double lat0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj = setRefPtProj(lon0, lat0);
    proj.type = CylEqDist;
    *projPtr = proj;
    return 1;
}

int GeogProjSetMercator(double lon0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj.type = Mercator;
    proj.params.lon0 = GeogLonR(lon0, 0.0);
    proj.rotation = 0.0;
    proj.cosr = 1.0;
    proj.sinr = 0.0;
    *projPtr = proj;
    return 1;
}

int GeogProjSetCylEqArea(double lon0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj.type = CylEqArea;
    proj.params.lon0 = GeogLonR(lon0, 0.0);
    proj.rotation = 0.0;
    proj.cosr = 1.0;
    proj.sinr = 0.0;
    *projPtr = proj;
    return 1;
}

int GeogProjSetLambertConfConic(double lon0, double lat0,
	double lat1, double lat2, struct GeogProj *proj_p)
{
    struct GeogProj proj;
    double n, F, rho0;
    double r0 = GeogREarth(NULL);

    proj.type = LambertConfConic;
    lat0 = GeogLatN(lat0);
    proj.params.LambertConfConic.lat0 = lat0;
    if ( proj.params.LambertConfConic.lat0 == 0.0 ) {
	return GeogProjSetMercator(lon0, proj_p);
    }
    proj.params.LambertConfConic.lon0 = GeogLonR(lon0, 0.0);
    if ( lat1 == lat2 ) {
	n = sin(lat1);
    } else {
	n = log(cos(lat1) / cos(lat2))
	    / log(tan(M_PI_4 + lat2 / 2.0) / tan(M_PI_4 + lat1 / 2.0));
    }
    if ( !isfinite(n) ) {
	fprintf(stderr, "  Lambert Conformal Conic n parameter not finite.\n");
	return 0;
    }
    F = cos(lat1) * pow(tan(M_PI_4 + lat1 / 2.0), n) / n;
    if ( !isfinite(F) ) {
	fprintf(stderr, "  Lambert Conformal Conic F parameter not finite.\n");
	return 0;
    }
    rho0 = r0 * F / pow(tan(M_PI_4 + lat0 / 2.0), n);
    proj.params.LambertConfConic.n = n;
    proj.params.LambertConfConic.F = F;
    proj.params.LambertConfConic.rho0 = rho0;
    proj.rotation = 0.0;
    proj.cosr = 1.0;
    proj.sinr = 0.0;
    *proj_p = proj;
    return 1;
}

int GeogProjSetLambertEqArea(double lon0, double lat0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj = setRefPtProj(lon0, lat0);
    proj.type = LambertEqArea;
    *projPtr = proj;
    return 1;
}

int GeogProjSetStereographic(double lon0, double lat0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj = setRefPtProj(lon0, lat0);
    proj.type = Stereographic;
    *projPtr = proj;
    return 1;
}

int GeogProjSetOrthographic(double lon0, double lat0, struct GeogProj *projPtr)
{
    struct GeogProj proj;

    proj = setRefPtProj(lon0, lat0);
    proj.type = Orthographic;
    *projPtr = proj;
    return 1;
}

/*
   Set members for a reference point based projection.
 */

static struct GeogProj setRefPtProj(double lon0, double lat0)
{
    struct GeogProj proj;

    lat0 = GeogLatN(lat0);
    proj.params.RefPt.lat0 = lat0;
    proj.params.RefPt.cos_lat0 = cos(lat0);
    proj.params.RefPt.sin_lat0 = sin(lat0);
    proj.params.RefPt.lon0 = GeogLonR(lon0, 0.0);
    proj.rotation = 0.0;
    proj.cosr = 1.0;
    proj.sinr = 0.0;
    return proj;
}

