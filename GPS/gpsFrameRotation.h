/*
 *  gpsFrameRotation.h
 *
 *  Created by Lee Jasper on Friday April 5, 2013.
 *  Copyright (c) 2013 Lee Jasper. All rights reserved.
 *
 *  Original MatLab software provided by Ben K. Bradley
 *  This is simply a conversion of Ben's code into FSW
 *   for the CICERO project.
 *  University of Colorado, Boulder
 *
 */

#ifndef _GPS_FRAME_ROTATION_H_
#define _GPS_FRAME_ROTATION_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "vector3D.h"
//#include "ciceroDefinitions.h"
#include "IAU2006_XYS.h"
#include "leapSecond.h"
//#include "RigidBodyKinematics.h"

//const float OMEGA_EARTH = 7.29211585530066e-5; //rad/s Rotation Rate of Earth

typedef struct {
    float    *array;
    int         used;
    int         size;
} fpArray;

typedef struct {
    float    **array;
    int         used;
    int         sizeX;
    int         sizeY;
} fpDoubleArray;

void wgs2gcrf(float_t *xEcef, float_t *vEcef, float_t *gpsTime, float_t *xEci, float_t *vEci);
void gcrf2wgs(float_t *xEci, float_t *vEci, float_t *gpsTime, float_t *xEcef, float_t *vEcef);
void gpstow2jd( float_t wn, float_t tow, float_t rollflag, float_t *jd_gps );
void format_JD( float_t *jd_orig, float_t *jd );
void jdtt2jdutc( float_t *jd_tt, float_t *jd_utc, int *tai_utc );
void ecef2eci_IAU2006CIOinterp( float_t *jd_tt, float_t *jd_utc, float_t *r0, float_t *r, float_t *v0, float_t *v, float_t IE[4][4] );
void eci2ecef_IAU2006CIOinterp(float_t *jd_tt, float_t *jd_utc, float_t *r0, float_t *r, float_t *v0, float_t *v, float_t EI[4][4]);
void computeERAmat( float_t *jd_ut1, float_t rERA[4][4], float_t *ERA );
void getXYs_simple( float_t mjd_tt, float_t *X, float_t *Y, float_t *s );
void interpLagrange( float_t xx, int p, int row0, float_t *yy);
void computeBPNmatrix( float_t X, float_t Y, float_t s, float_t BPN[4][4]);

void gcrf2j2000(float_t *xEci, float_t *vEci, float_t *xJ2000, float_t *vJ2000, float_t B[4][4]);
void j20002gcrf(float_t *xJ2000, float_t *vJ2000, float_t *xEci, float_t *vEci, float_t B[4][4]);

void initArray(fpArray *a, int initialSize);
void insertArray(fpArray *a, float_t element);
void freeArray(fpArray *a);
void initDoubleArray(fpDoubleArray *a, int initialDim1, int initialDim2);
void Mult1xN_NxM( fpArray *a, fpDoubleArray *b, float_t *c );
void freeDoubleArray(fpDoubleArray *a, int dim1);



#endif
