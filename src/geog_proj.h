/*
   -	geog_proj.h --
   -
   -		This header file defines constants and declares structures and
   -		functions that manage and apply cartographic projections.
   .
   .	Copyright (c) 2012, Gordon D. Carrie. All rights reserved.
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
   .	$Revision: 1.3 $ $Date: 2014/01/24 22:43:39 $
 */

/*

   Ref.
   Snyder, John P.
   Map Projections used by the U.S. Geological Survey.
   (Geological Survey bulletin ; 1532)
   United States Government Printing Office, Washington:  1982.

 */

#ifndef _GEOGPROJ_H_
#define _GEOGPROJ_H_

/* Identifiers for the currently recognized projection types */ 
enum ProjType {
    CylEqDist, CylEqArea, Mercator, LambertConfConic,
    LambertEqArea, Orthographic, Stereographic
};

/*
   Structures of this store parameters and functions that making conversions
   between geographic and map coordinates. Clients should not refer to fields
   in this structure directly.
 */

struct GeogProj {
    enum ProjType type;			/* Projection type */
    union {
	double lon0;			/* Longitude of map origin */
	struct {
	    double lon0, lat0;		/* Longitude, latitude or map origin */
	    double cos_lat0;		/* Cosine of lat0 */
	    double sin_lat0;		/* Sine of lat0 */
	} RefPt;
	struct {
	    double lon0;		/* Longitude of map origin, vertical
					   meridian if rotation == 0.0 */
	    double lat0;		/* Latitude of map origin */
	    double lat1, lat2;		/* Standard parallels */
	    double rho0;		/* See Snyder, p. 105 */
	    double n;			/* See Snyder, p. 105 */
	    double F;			/* See Snyder, p. 105 */
	} LambertConfConic;
    } params;
    double rotation;			/* Rotation angle (clockwise degrees).
					   This specifies the angle between
					   geographic north and north on the
					   map. It can also be regarded as the
					   orientation of the surface onto which
					   the Earth's surface is projected */
    double cosr, sinr;			/* Cosine and sine of rotation */
};

int GeogProjXYToLonLat(double, double, double *, double *, struct GeogProj *);
int GeogProjLonLatToXY(double, double, double *, double *, struct GeogProj *);
int GeogProjSetCylEqDist(double, double, struct GeogProj *);
int GeogProjSetCylEqArea(double, struct GeogProj *);
int GeogProjSetMercator(double, struct GeogProj *);
int GeogProjSetLambertConfConic(double, double, double, double,
	struct GeogProj *);
int GeogProjSetLambertEqArea(double, double, struct GeogProj *);
int GeogProjSetStereographic(double, double, struct GeogProj *);
int GeogProjSetOrthographic(double, double, struct GeogProj *);
void GeogProjSetRotation(struct GeogProj *, double);
int GeogProjSetFmStr(char *, struct GeogProj *);

#endif
