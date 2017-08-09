/*
 *  gpsFrameRotation.c
 *
 *  Created by Lee Jasper on Friday June 13, 2013.
 *  Copyright (c) 2013 Lee Jasper. All rights reserved.
 *
 *  Original MatLab software provided by Ben K. Bradley.
 *  This is simply a conversion of Ben's code and the FSW
 *   from7 the CICERO project into the FSW for the ALL-STAR
 *	 program.
 *  University of Colorado, Boulder
 *
 */

#include "gpsFrameRotation.h"
const float OMEGA_EARTH = 7.29211585530066e-5;

/*
 * Function: wgs2gcrf
 * Purpose: Computes the inertial position of a satellite (GCRF) given GPS data
 * Author: Ben K. Bradley
 * Date: 03/22/2013
 * Inputs:
 *  xEcef       = Position vector, km, in GPS (ITRF, ECEF)
 *  vEcef       = Velocity vector, km/s, in GPS (ITRF, ECEF)
 *  gpsTime[2]  = GPS week number and seconds into GPS week (from 0hr on Sunday)
 * Outputs:
 *  xEci        = Position vector, km, in GCRF (ECI)
 *  vEci        = Velocity vector, km/s, in GCRF (ECI)
 */
void wgs2gcrf(float_t *xEcef, float_t *vEcef, float_t *gpsTime, float_t *xEci, float_t *vEci)
{
    int         rollFlag = 1;
    int         tai_utc;
    float_t    jd_gps[2];
    float_t    jd_tt[2];
    float_t    jd_utc[2];
    float_t    IE[4][4];

    /* Convert GPS time of week to Julian Date */
    gpstow2jd(gpsTime[0], gpsTime[1], rollFlag, jd_gps);

    /* Terrestrial Time (TT) Julian Date */
    jd_tt[0] = jd_gps[0];
    jd_tt[1] = jd_gps[1];
    /* 19 is number of leap seconds in 1980, start of GPS */
    /* 32.184 is conversion from TAI to TT */
    jd_tt[1] = (19.0 + 32.184)/86400.0 + jd_tt[1];

    /* UTC Julian date */
    jdtt2jdutc(jd_tt, jd_utc, &tai_utc);

    /* Error checking: jd_utc should not be 0 */
    if(jd_utc[0] == 0 && jd_utc[1] == 0) {
        xEci[1] = 0;
        xEci[2] = 0;
        xEci[3] = 0;
        return;
    }

    /* Frame Transformation */
    /* Convert ECEF to ECI */
    ecef2eci_IAU2006CIOinterp(jd_tt, jd_utc, xEcef, xEci, vEcef, vEci, IE);
}



/*
 * Function: gcrf2wgs
 * Purpose: Computes the GPS position given inertial (GCRF) position
 * Author: Ben K. Bradley
 * Date: 03/22/2013
 * Inputs:
 *  xEci        = Position vector, km, in GCRF (ECI)
 *  vEci        = Velocity vector, km/s, in GCRF (ECI)
 *  gpsTime[2]  = GPS week number and seconds into GPS week (from 0hr on Sunday)
 * Outputs:
 *  xEcef       = Position vector, km, in GPS (ITRF, ECEF)
 *  vEcef       = Velocity vector, km/s, in GPS (ITRF, ECEF)
 */
void gcrf2wgs(float_t *xEci, float_t *vEci, float_t *gpsTime, float_t *xEcef, float_t *vEcef){
    int         rollFlag;
    int         tai_utc;
    float_t    jd_gps[2];
    float_t    jd_tt[2];
    float_t    jd_utc[2];
    float_t    IE[4][4];

    rollFlag = 1;
    
    /* Convert GPS time of week to Julian Date */
    gpstow2jd( gpsTime[0], gpsTime[1], rollFlag, jd_gps );

    /* Terrestrial Time (TT) Julian Date */
    jd_tt[0] = jd_gps[0];
    jd_tt[1] = jd_gps[1];
    /* 19 is number of leap seconds in 1980, start of GPS */
    /* 32.184 is conversion from TAI to TT */
    jd_tt[1] = (19.0 + 32.184)/86400.0 + jd_tt[1];

    /* UTC Julian date */
    jdtt2jdutc(jd_tt, jd_utc, &tai_utc);

    /* Error checking: jd_utc should not be 0 */
    if( jd_utc[0] == 0 && jd_utc[1] == 0 ){
        xEci[1] = 0;
        xEci[2] = 0;
        xEci[3] = 0;
        return;
    }

    /* Frame Transformation */
    /* Convert ECEF to ECI */
    eci2ecef_IAU2006CIOinterp( jd_tt, jd_utc, xEci, xEcef, vEci, vEcef, IE );
}




/*
 * Function: gpstow2jd
 * Purpose: Calculates a Julian date from a GPS week number and time of week
 * Author: Ben K. Bradley
 * Date: 04/28/2010
 * Inputs:
 *  wn          = GPS week number
 *  tow         = seconds into GPS week (starting from 0hr on Sunday)
 *  rollflag    = flag dictating epoch for wn
 *                  1 = 1980.01.06
 *                  2 = 1999.08.22
 * Outputs:
 *  jd_gps[2]   = 2-part julian date [day; dayfraction]
 */
void gpstow2jd(float_t wn, float_t tow, float_t rollflag, float_t *jd_gps)
{
    float_t jd_gps_epoch;
    float_t jd_gps_temp[2];

    /* Set both GPS week number origins */
    /* 06jan80 = 2444244.5, Julian Day = January 6, 1980 */
    /* 22aug99 = 2451412.5, Julian Day = August 22, 1999 */
    if( rollflag == 1 ){
        jd_gps_epoch = 2444244.5;
    }
    else if ( rollflag == 2 ) {
        jd_gps_epoch = 2451412.5;
    }
    else {
        /* Error occurred, this will trigger bad data */
        jd_gps_epoch = 0.0; 
    }

    /* Compute Julian Date */
    jd_gps_temp[0] = wn*7.0 + jd_gps_epoch;
    jd_gps_temp[1] = tow/86400.0;

    /* Enforce correct [day dayfraction] formatting */
    format_JD(jd_gps_temp, jd_gps);
}

/*
 * Function: format_JD
 * Purpose: This takes a 2-part Julian Date and makes sure that it is in the 
 *  correct day and day fraction format. This format is required for cases 
 *  requiring high precision ( e.g. JD = [2451522.5  0.2134] ).  This will also
 *  work when the input is full Julian Date in the first cell and a zero in the 
 *  second cell ( e.g. JD = [2456161.1158  0] ).  Additionally, this algorithm 
 *  will work when the second cell contains multiple days 
 *  (e.g. JD = [2451522.5 4.2234]).
 * Author: Ben K. Bradley
 * Date: 08/21/2012
 * Inputs:
 *  jd_orig[2]  = 2-part Julian Date [day; dayfraction]
 * Outputs:
 *  jd          = 2-part Julian Date [day; dayfraction]
 */
void format_JD(float_t *jd_orig, float_t *jd)
{
    float_t jd1;
    float_t extra2;
    float_t add_days;

    jd[0] = jd_orig[0];
    jd[1] = jd_orig[1];

    /* Store extra time in and modify JD */
    jd1    = jd_orig[0] - 0.5;
    extra2 = jd1 - floor(jd1);

    jd[0] = floor(jd1) + 0.5;
    jd[1] = jd[1] + extra2;

    /* Make sure second array value is between 0 and 1 */
    if((jd[1] >= 1) || (jd[1] < 0)) {
        add_days = floor(jd[1]);
        jd[0] = jd[0] + add_days;
        jd[1] = jd[1] - add_days;
    }
}

/*
 * Function: jdtt2jdutc
 * Purpose: Calculates a UTC Julian Date from a TT (Terrestrial Time) Julian 
 *  Date. This function uses the leapsec.dat file to determine the number of 
 *  leap seconds to apply. Please note that when the resulting UTC time is on 
 *  the day a leap second is added, the Julian Date is not 100% correct.
 *  However, the resulting Julian Date will be perfectly reasonable for use
 *  in interpolating EOP data, for instance.
 * Author: Ben K. Bradley
 * Date: 08/22/2012
 * Inputs:
 *  jd_tt[2]   = 2-part Julian Date (TT) [day; dayfraction]
 * Outputs:
 *  jd_utc[2]  = 2-part Julian Date (UTC) [day; dayfraction]
 *  tai_utc    = Number of leap seconds used to convert TAI to UTC
 */
void jdtt2jdutc( float_t *jd_tt, float_t *jd_utc, int *tai_utc )
{
    float_t		jd_tai[2];
    float_t		mjd_tai;
    float_t		mjd_utc;
    int			sizeLeap;
    int         row;
    int         i;
    int         leap;
    int         numCols;

    /* Leap second has 6 'columns' per leap second update. This is actually 
     * stored as a single array so offsets are needed */
    numCols = 6;

    /* Compute International Atomic Time (TAI) */
    jd_tai[0] = jd_tt[0];
    jd_tai[1] = jd_tt[1];
    jd_tai[1] = -32.184/86400.0 + jd_tt[1];
    mjd_tai = jd_tai[1] + ( jd_tai[0] - 2400000.5 ) ;
    
    /* TAI date/time is before official start of TAI time: Jan 1 1972, 0hr */
    if( mjd_tai < 41317 ) {
        jd_utc[0] = 0;
        jd_utc[1] = 0;
        *tai_utc = 0;
        return;
    }

    /* Use leap seconds */
    sizeLeap = sizeof(leapsecond)/sizeof(*leapsecond);
    /* Find number of leap seconds */
    /* This leap sec count is best guess given TAI time */
    row = 0;
    for(i = 0; i < sizeLeap; i = i + numCols) {
        if(mjd_tai >= leapsecond[i + 3]) {
            row = i;
        }
    }

    if(row == 0) {
        leap = 0;
    } else {
        /* TAI-UTC */
        leap = (int)leapsecond[row + 5];
    }

    /* Convert to UTC */
    jd_utc[0] = jd_tai[0];
    jd_utc[1] = jd_tai[1];

    jd_utc[1] = (float_t)-leap/86400.0 + jd_tai[1];
    mjd_utc = jd_utc[1] + ( jd_utc[0] - 2400000.5 );

    /*Check if current UTC guess crossed a leap second boundary */
    if( mjd_utc < leapsecond[row + 3] ){
        /* Boundary was crossed. Leap second count is one less */
        leap = leap - 1;
        jd_utc[1] = (float_t)-leap/86400 + jd_tai[1];

        format_JD( jd_utc,  jd_utc);

        /* If this boundary crossing and leap second modification results in 
        the time being between 0 <= t < 1 seconds into the day, this means 
        that the correct UTC time would actually have fallen on the 60th 
        second if we were representing this time as a date/time vector. Since 
        the Julian date cannot accurately represent the 60th second, we at 
        least force this resulting UTC Julian day to fall just before the leap
        second is introduced.  This will allow for correct interpolation of
        EOPs (specifically dUT1) when using this UTC Julian date. */

        if((jd_utc[1] >= (float_t)0.0) && (jd_utc[1] < (float_t)1/86400)) {
            jd_utc[1] = (86400.0 + jd_utc[1]*86400.0)/86401.0;
        }
    }

    /* Make sure JD is in proper 2-part format */
    format_JD(jd_utc, jd_utc);

    *tai_utc = leap;
}

/*
 * Function: ecef2eci_IAU2006CIOinterp
 * Purpose: Converts a position vector in the ITRF(ECEF) frame to the GCRF(ECI)
 *  frame using the IAU 2006 precession and IAU 2000A_R06 nutation theories. 
 *  This routine employs a hybrid of the "Full Theory" using Fukushima-Williams
 *  angles and the CIO-based method.  
 *
 *  Specifically, this routine interpolates a table of X,Y,s values and then
 *  uses them to construct the BPN matrix directly.  The X,Y,s values in the
 *  data table were generated using Fukushima-Williams angles and the 
 *  IAU 2000A_R06 nutation theory.  This general scheme is outlined in [3]
 *  and [4]. Earth orientation parameters (EOPs) are not used for speed
 *  and simplicity purposes. Increased accuracy can be gained if the UT1
 *  Julian date is passed in instead of UTC as specified.  
 *
 *  The function init_XYs2006.m must be run before this function is called.  
 *  It loads data files necessary for the transformation. 
 * Author: Ben K. Bradley
 * Date: 03/22/2013
 * Inputs:
 *  jd_tt[2]    = 2-part Julian Date (TT) [day; dayfraction]
 *  jd_utc[2]   = 2-part Julian Date (UTC)[day; dayfraction]
 *  r0          = Position vector in ITRF (ECEF)
 *  v0          = Velocity vector in ITRF (ECEF)
 * Outputs:
 *  r           = position vector in GCRF (ECI)
 *  v           = velocity vector in GCRF (ECI)
 *  IE          = rotation matrix used to convert position from ecef to eci
 */
void ecef2eci_IAU2006CIOinterp(float_t *jd_tt, float_t *jd_utc, 
                               float_t *r0, float_t *r, 
                               float_t *v0, float_t *v, 
                               float_t IE[4][4])
{
    float_t		R[4][4];
    float_t		BPN[4][4];
    float_t		X;
    float_t		Y;
    float_t		s;
    float_t		ERA;
    float_t		v_E[4];
    float_t		wXr[4];
    float_t		omegaEI[4];

    /* Construct Earth rotation angle matrix, R (TIRS to CIRS) */
    /* e.g. r_cirs = R * r_tirs */
    computeERAmat(jd_utc, R, &ERA);

    /*Construct bias-precession-nutation matrix BPN (CIRS->GCRS(ICRS)) */
    getXYs_simple(jd_tt[1] + (jd_tt[0] - 2400000.5), &X, &Y, &s);

    /* Error Checking: correct dates */
    if( X == 0 && Y == 0 && s == 0 ) {
        set3(0,0,0,r);
        return;
    }

    /* Calculate BPN matrix */
    computeBPNmatrix( X, Y, s, BPN);

    /* Transform Position Vector */
    /* FI = BPN * R; */
    MdotM(BPN, R, IE);
    /* r = IE * r0; */
    Mdot(IE, r0, r);

    /* Transform Velocity Vector */
    /* v_E = E d/dt(r0) + omega_E/I x r0 = v0 + omega_E/I x r0 */
    set3(0, 0, OMEGA_EARTH, omegaEI);
    cross(omegaEI, r0, wXr);
    add(v0, wXr, v_E);

    /* v_I = [IE] v_E */
    Mdot(IE, v_E, v);
}



/*
 * Function: eci2ecef_IAU2006CIOinterp
 * Purpose: Converts a position vector in the ITRF(ECEF) frame to the GCRF(ECI)
 *  frame using the IAU 2006 precession and IAU 2000A_R06 nutation theories. 
 *  This routine employs a hybrid of the "Full Theory" using Fukushima-Williams 
 *  angles and the CIO-based method.  
 *
 *  Specifically, this routine interpolates a table of X,Y,s values and then
 *  uses them to construct the BPN matrix directly.  The X,Y,s values in the
 *  data table were generated using Fukushima-Williams angles and the 
 *  IAU 2000A_R06 nutation theory.  This general scheme is outlined in [3]
 *  and [4]. Earth orientation parameters (EOPs) are not used for speed
 *  and simplicity purposes. Increased accuracy can be gained if the UT1
 *  Julian date is passed in instead of UTC as specified.  
 *
 *  The function init_XYs2006.m must be run before this function is called.  
 *  It loads data files necessary for the transformation.
 * Author: Ben K. Bradley
 * Date: 03/22/2013
 * Inputs:
 * Inputs:
 *  jd_tt[2]    = 2-part Julian Date (TT) [day; dayfraction]
 *  jd_utc[2]   = 2-part Julian Date (UTC)[day; dayfraction]
 *  r0          = Position vector in GCRF (ECI)
 *  v0          = Velocity vector in GCRF (ECI)
 * Outputs:
 *  r           = position vector in ITRF (ECEF)
 *  v           = velocity vector in ITRF (ECEF)
 *  IE          = rotation matrix used to convert position from ecef to eci
 */
void eci2ecef_IAU2006CIOinterp(float_t *jd_tt, float_t *jd_utc, 
                               float_t *r0, float_t *r, 
                               float_t *v0, float_t *v, 
                               float_t EI[4][4])
{
    float_t R[4][4];
    float_t BPN[4][4];
    float_t IE[4][4];
    float_t X;
    float_t Y;
    float_t s;
    float_t ERA;
    float_t v_E[4];
    float_t wXr[4];
    float_t omegaEI[4];
    
    /* Construct Earth rotation angle matrix, R (TIRS to CIRS) */
    /* e.g. r_cirs = R * r_tirs */
    computeERAmat( jd_utc, R, &ERA );

    /*Construct bias-precession-nutation matrix BPN (CIRS->GCRS(ICRS)) */
    getXYs_simple( jd_tt[1] + (jd_tt[0] - 2400000.5), &X, &Y, &s );

    /* Error Checking: correct dates */
    if( X == 0 && Y == 0 && s == 0 ) {
        set3(0,0,0,r);
        return;
    }

    /* Calculate BPN matrix */
    computeBPNmatrix( X, Y, s, BPN);

    /* Transform Position Vector */
    /* IE = BPN * R; */
    MdotM(BPN, R, IE);
    transpose(IE, EI);
    /* r = [EI] * r0; */
    Mdot(EI, r0, r);
    
    /* Transform Velocity Vector */
    /* v_E = [EI] v_I */
    Mdot(EI, v0, v_E);

    /* E d/dt(r) = v_E - omega_E/I x r */
    set3(0, 0, 7.292115146706979e-5, omegaEI);
    cross(omegaEI, r, wXr);
    sub(v_E, wXr, v);
}



/*
 * Function: compute_ERAmat
 * Purpose: Computes the Earth Rotation Angle (ERA) and the ERA rotation matrix
 *  based on UT1 time. The ERA is modulated to lie within [0,2*pi] and is 
 *  computed using the precise equation given by Eq. 5.15 in [1].
 *
 *  The ERA is the angle between the Celestial Intermediate Origin, CIO, and 
 *  Terrestrial Intermediate Origin, TIO (a reference meridian 100m offset 
 *  from Greenwich meridian).
 * Author: Ben K. Bradley
 * Date: 05/14/2012
 * Inputs:
 *  jd_tt[2]    = 2-part Julian Date (UT1) [day; dayfraction]
 * Outputs:
 *  rERA        = ERA rotation matrix
 *  ERA         = Earth Rotation Angle (0 <= ERA <= 2*pi)
 */
void computeERAmat(float_t *jd_ut1, float_t rERA[4][4], float_t *ERA)
{
    float_t		F;
    float_t		era;
    float_t		ce;
    float_t		se;
    float_t		twoPi = 2 * M_PI;

    /* Compute ERA */
    /* F = (jd_ut1[0] % 1) + (jd_ut1[1] % 1); */
    F = (jd_ut1[0] - floor(jd_ut1[0])) + (jd_ut1[1] - floor(jd_ut1[1]));

    era = (float_t)(twoPi * (F + 0.7790572732640 + 0.00273781191135448 * 
        (jd_ut1[1] + (jd_ut1[0] - 2451545))));

    /* era = era % twoPi; */
	/* TODO: fmodf MAY NOT EXIST. INCLUDE MATH.H ? */
    era = fmodf(era, twoPi);

    if( era < 0 ) {
        era = era + twoPi;
    }
    
    /* Construct ERA rotation matrix */
    /* R = rot3(-ERA) */
    ce = cos(era);
    se = sin(era);
    setMatrix(ce, -se, 0,
              se, ce, 0,
              0, 0, 1, rERA);
    *ERA = era;
}


/*
 * Function: getXYs_simple
 * Purpose: 
 *  Interpolates X,Y, and s loaded by init_XYs2006.m using an 11th-order 
 *  Lagrange interpolation method. The init_XYsdata.m function must be 
 *  called before get_XYs is used.  This function uses the XYs data set that 
 *  has been loaded as a matrix. Each of the three values listed below are 
 *  tabulated at 0h TT of each day.
 *  
 *  X: x-coordinate of the Celestial Intermediate Pole (CIP)
 *  Y: y-coordinate of the Celestial Intermediate Pole (CIP)
 *  s: Celestial Intermediate Origin (CIO) locator
 * Author: Ben K. Bradley
 * Date: 05/14/2012
 * Inputs:
 *  mjd_tt  = Modified Julian Date (TT)
 * Outputs:
 *  X       = X-coord of the CIP in the GCRF (part of unit vector) (radians)
 *  Y       = Y-coord of the CIP in the GCRF (part of unit vector) (radians)
 *  s       = CIO Locator (radians)
 *   
 */
void getXYs_simple(float_t mjd_tt, float_t *X, float_t *Y, float_t *s)
{
    float_t		sec2rad;
    float_t		tempans[3];
    int         seed;
    int         numCols;
    int         numElements;

    /* indicies in array that should be interpreted as columns */
    numCols = 7;
    numElements = sizeof(iau2006_XYs)/sizeof(*iau2006_XYs);

    sec2rad = 4.848136811095359935899141e-6; // arcseconds to radians

    /* Interpolate XYs data */

    /*Date is within data range: checking Julian date */
    if( (mjd_tt >= iau2006_XYs[0 + 3]) & (mjd_tt <= iau2006_XYs[numElements-4]) ) {
        /* Eliminates need for locate_bisection call. We can do this because
            we know iau2006_XYs is tabulated at 1 day steps */
        seed = (int)(floor(mjd_tt) - iau2006_XYs[0 + 3] + 1);

        /*Interpolate using 11th-order Lagrange */
        interpLagrange(mjd_tt, 7, seed, tempans);

        *X = (float_t)(tempans[0] * sec2rad);
        *Y = (float_t)(tempans[1] * sec2rad);
        *s = (float_t)(tempans[2] * sec2rad);
    } else {
        /*Date is outside of XYs table. This should not happen */
        *X = 0;
        *Y = 0;
        *s = 0;
    }
}

/*
 * Function: interpLagrange
 * Purpose: Interpolates data using the Langrange method of order p
 * Author: Ben K. Bradley
 * Date: 05/25/2011 (updated: 8/28/2012)
 * Inputs:
 *  xx      = single x value to interpolate at
 *  p       = order of interpolation (e.g. 4,5,7,9,...)
 *  row0    = the row number of X that locate_bisection.m
 *            would return. i.e. row of X such that   
 *            X(row0) <= xx < X(row0+1)
 *            Input a -1 if you don't know. locate_bisection will
 *            be called then to determine the correct reference row
 * Outputs:
 *  yy      = interpolated Y value(s)
 */
void interpLagrange( float_t xx, int p, int row0, float_t *yy)
{
    int             N;
    int             numElements;
    int             numCols;
    int             iStart;
    int             iEnd;
    int             i;
    int             j;
    float_t			No2;
    int             nn;
    int             numRows;
    float_t			xrow0;
    float_t			pjValue;
    fpArray         X;
    fpArray         Pj;
    fpDoubleArray   Y;

    /* indicies in array that should be interpreted as columns */
    numCols = 7;
    numElements = sizeof(iau2006_XYs)/sizeof(*iau2006_XYs);
    numRows = numElements/numCols;
    
    /* Number of data points to use for interpolation (e.g. 8, 9, 10,...) */
    N = p + 1;

    /* Initialize arrays */
    initArray(&X, N);
    initArray(&Pj, N);
    initDoubleArray(&Y, N, 3); /*Always cols 5,6,7 */

    if(numRows < N) {
        /*Not enough data points for Lagrange Interpolation */
        yy[0] = 0;
        yy[1] = 0;
        yy[2] = 0;
        return;
    }

    /* Compute number of elements on either side of middle element to grab */
	/* TODO: need floor defined */
    No2 = 0.5 * N;
    nn = (int)floor(No2);

    /* Unnecessary in flight code */
    //if( row0 == -1) {
    //    /*Find row of indep. variable such that X(row0) <= xx < X(row0+1) */
    //    row0 = locateBisection(xx);
    //    /* Note: row0 = 0 when xx is before all data in X
    //             row0 = length(X) when xx is after all data in X */
    //}

    if(No2 - nn == 0) {
        /*adjust row0 in case near data set endpoints*/
        if( (N == numRows) || (row0 < No2) ) {
            row0 = (int)No2;
        } 
        else if ( row0 > (numRows - No2) ) {
            row0 = (int)(numRows - No2);
        }

        /* Trim to relevant data points */
        iStart = (int)(row0 - No2)*numCols;
        iEnd = (int)(row0 + No2)*numCols;
        j = 0;
        for(i = iStart; i < iEnd; i=i+numCols) {
            insertArray( &X, iau2006_XYs[i+3] );
            Y.array[j][0] = iau2006_XYs[i+4];
            Y.array[j][1] = iau2006_XYs[i+5];
            Y.array[j][2] = iau2006_XYs[i+6];
            j++;
        }
    } else {
        if((N == numRows) || (row0 < nn+1)) {
            row0 = nn + 1;
        } else if(row0 == numRows) {
            row0 = numRows - p;
        } else {
            xrow0 = iau2006_XYs[(row0-1)*numCols + 3]; 
            if ( xx - xrow0 > 0.5 ) {
                row0 = row0 + 1;
            }
        }
        /* Trim to relevant data points */
        iStart = (row0 - nn - 1)*numCols;
        iEnd = (row0 + nn)*numCols;
        for(i = iStart; i < iEnd; i++){
            insertArray( &X, iau2006_XYs[i+3] );
            Y.array[i][0] = iau2006_XYs[i+4];
            Y.array[i][1] = iau2006_XYs[i+5];
            Y.array[i][2] = iau2006_XYs[i+6];
        }
    }

    /* Loop over all data points being included */
    for(j = 0; j < N; j++) {
        insertArray( &Pj, 1 );
        for(i = 0; i < N; i++) {
            if( j != i ) {
                pjValue = Pj.array[j] * (-xx + X.array[i])/(-X.array[j] + 
                    X.array[i] );
                Pj.array[j] = pjValue;
            }
        }
    }

    /* yy = Pj * Y */
    Mult1xN_NxM( &Pj, &Y, yy );

    freeArray(&X);
    freeArray(&Pj);
    freeDoubleArray(&Y, N);
}

/*
 * Function: computeBPNmatrix
 * Purpose: Computes the Bias-Precession-Nutation matrix required for the 
 *  CIO-based transformation between the GCRF/ITRF frames. The Z-coordinate
 *  of the CIP is also computed and available as an output. This is a slightly 
 *  simplified routine that assumes cos(s)=1 and sin(s)=s.
 * Author: Ben K. Bradley
 * Date: 05/20/2012
 * Inputs:
 *  X       = X-coord of the CIP in the GCRF (part of unit vector) (radians)
 *  Y       = Y-coord of the CIP in the GCRF (part of unit vector) (radians)
 *  s       = CIO Locator (radians)
 * Outputs:
 *  BPN     = Bias-Precession-Nutation matrix, e.g. r_gcrs = BPN * r_cirs
 *   
 */
void computeBPNmatrix( float_t X, float_t Y, float_t s, float_t BPN[4][4])
{
    float_t    aa;
    float_t    fXY[4][4];
    float_t    fs[4][4];

    aa = 1/(1 + sqrt(1 - X*X - Y*Y) );

    setMatrix(1-aa*X*X, -aa*X*Y, X,
                -aa*X*Y, 1-aa*Y*Y, Y,
                -X, -Y, 1-aa*(X*X + Y*Y), fXY);

    /* Approximations: cos(s) = s and sin(s) = 1 */
    setMatrix(1, s, 0,
              -s, 1, 0,
              0, 0, 1, fs);

    MdotM(fXY, fs, BPN);

}

/*
 * Function: gcrf2j2000
 * Purpose: Computes the Frame Bias matrix 'B' between GCRF and J2000. The 
 *  rotation is a constant rotation and is very small, on the order of 4 
 *  meters difference. The values have been obtained from the Astronomical 
 *  Almanac 2006, sect. B27
 * Author: Lee E. Z. Jasper
 * Date: 05/23/2013
 * Inputs:
 *  xEci    = position in GCRF reference frame
 *  vEci    = velocity in GCRF reference frame
 * Outputs:
 *  xJ2000  = position in J2000 reference frame
 *  vJ2000  = velocity in J2000 reference frame
 *  B       = Frame Bias matrix, e.g. r_gcrf = B * r_j2000
 */
void gcrf2j2000(float_t *xEci, float_t *vEci, 
                float_t *xJ2000, float_t *vJ2000, float_t B[4][4])
{
    /* 
    // Computations for matrix. Since matrix is constant it is simply inputted: 
    float_t    arcsec2rad;
    float_t    delta_a0;
    float_t    delta_psiB;
    float_t    epsilon0;
    float_t    eta0;
    float_t    zeta0;
    float_t    r3[4][4];
    float_t    r2[4][4];
    float_t    r1[4][4];
    float_t    temp[4][4];

    arcsec2rad = 4.8481368e-6; //[arcsec]->[rad]

    // Constants from Astro. Almanic 2006
    delta_a0 = -0.0146 * arcsec2rad;
    delta_psiB = -0.041775 * arcsec2rad;
    epsilon0 = 84381.448 * arcsec2rad;
    eta0 = -0.0068192 * arcsec2rad;

    zeta0 = delta_psiB*sin(epsilon0);

    // Frame Bias (from J2000 to GCRF)
    //    B = EA_DCM(-delta_a0,3)*EA_DCM(-zeta0,2)*EA_DCM(eta0,1); 
    Euler3(-delta_a0, r3);
    Euler2(-zeta0, r2);
    Euler1(eta0, r1);

    MdotM(r1, r2, temp);
    MdotM(temp, r3, B_trans);
    */

    /* From GCRF -> J2000 */
    setMatrix(0.9999999999999942, -0.0000000707827974, 0.0000000805621715,
     0.0000000707827948, 0.9999999999999969, 0.0000000330604145,
    -0.0000000805621738, -0.0000000330604088, 0.9999999999999962, B);

    Mdot(B, xEci, xJ2000);
    Mdot(B, vEci, vJ2000);
}

/*
 * Function: j20002gcrf
 * Purpose: Computes the Frame Bias matrix 'B' between GCRF and J2000. The 
 *  rotation is a constant rotation and is very small, on the order of 4 
 *  meters difference. The values have been obtained from the Astronomical 
 *  Almanac 2006, sect. B27
 * Author: Lee E. Z. Jasper
 * Date: 05/23/2013
 * Inputs:
 * Inputs:
 *  xJ2000  = position in J2000 reference frame
 *  vJ2000  = velocity in J2000 reference frame
 * Outputs:
 *  xEci    = position in GCRF reference frame
 *  vEci    = velocity in GCRF reference frame
 *  B       = Frame Bias matrix, e.g. r_gcrf = B * r_j2000
 */
void j20002gcrf(float_t *xJ2000, float_t *vJ2000, 
                float_t *xEci, float_t *vEci, float_t B[4][4])
{
    float_t    B_trans[4][4];
    /* 
    // Computations for matrix. Since matrix is constant it is simply inputted: 
    float_t    arcsec2rad;
    float_t    delta_a0;
    float_t    delta_psiB;
    float_t    epsilon0;
    float_t    eta0;
    float_t    zeta0;
    float_t    r3[4][4];
    float_t    r2[4][4];
    float_t    r1[4][4];
    float_t    temp[4][4];

    arcsec2rad = 4.8481368e-6; //[arcsec]->[rad]

    // Constants from Astro. Almanic 2006
    delta_a0 = -0.0146 * arcsec2rad;
    delta_psiB = -0.041775 * arcsec2rad;
    epsilon0 = 84381.448 * arcsec2rad;
    eta0 = -0.0068192 * arcsec2rad;

    zeta0 = delta_psiB*sin(epsilon0);

    // Frame Bias (from J2000 to GCRF)
    //    B = EA_DCM(-delta_a0,3)*EA_DCM(-zeta0,2)*EA_DCM(eta0,1); 
    Euler3(-delta_a0, r3);
    Euler2(-zeta0, r2);
    Euler1(eta0, r1);

    MdotM(r1, r2, temp);
    MdotM(temp, r3, B_trans);
    */

    /* From GCRF -> J2000 */
    setMatrix(0.9999999999999942, -0.0000000707827974, 0.0000000805621715,
     0.0000000707827948, 0.9999999999999969, 0.0000000330604145,
    -0.0000000805621738, -0.0000000330604088, 0.9999999999999962, B_trans);

    /*From J2000 -> GCRF */
    transpose(B_trans, B);

    Mdot(B, xJ2000, xEci);
    Mdot(B, vJ2000, vEci);
}

void initArray(fpArray *a, int initialSize) {
  a->array = (float_t *)malloc(initialSize * sizeof(float_t));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(fpArray *a, float_t element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (float_t *)realloc(a->array, a->size * sizeof(float_t));
  }
  a->array[a->used++] = element;
}

void freeArray(fpArray *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void initDoubleArray(fpDoubleArray *a, int initialDim1, int initialDim2) {
    int i;

    a->array = (float_t**)malloc(initialDim1 * sizeof(float_t*));
    for (i = 0; i < initialDim1; i++) {
        a->array[i] = (float_t *)malloc(initialDim2 * sizeof(float_t));
    }
    a->used = 0;
    a->sizeX = initialDim1;
    a->sizeY = initialDim2;
}

void Mult1xN_NxM( fpArray *a, fpDoubleArray *b, float_t *c ) {
    int i;
    int j;
    int N;
    int N1;
    int M;

    N = a->size;
    N1 = b->sizeX;
    M = b->sizeY;

    if( N != N1 ) {
        /* Dimension mis-match */
        return;
    }

    for( i = 0; i < M; i++ ){
        c[i] = 0;
        for( j = 0; j < N; j++ ){
            c[i] = c[i] + (a->array[j] * b->array[j][i] );
        }
    }
}

void freeDoubleArray(fpDoubleArray *a, int dim1) {
    int i;

    for (i = 0; i < dim1; i++) {
        free(a->array[i]);
    }

    free(a->array);
    a->array = NULL;
    a->used = a->sizeX = a->sizeY = 0;
}
