/*
 * GPS.h
 *
 * Created: 7/11/2013 4:19:40 PM
 *  Author: Ryan Montoya
 * Contains relevant data to GPS: ECEF and ECIF conversions in a struct
 */ 

#ifndef GPS_H_
#define GPS_H_

#include <stdint.h>
/* This includes the previous CDH protocols, might need to be updated once CDH protocols are final.
#include "CDH/BufferPool.h"
*/

//GPS struct
typedef struct
{
	/* GPS Time */
	uint32_t GPS_Week;		// GPS Week number
	float GPS_TOW;			// Real GPS Time Of Week [s] obtained from the GPS unit. It is update automatically through Opcodes.
	float last_GPS_TOW;		// Stores the previous TOW so we can compare and determine if there is a new GPS solution
	float sysTime;			// Time since the last GPS solution [s]. Reset to zero every time we get a new GPS solution.
	int   offset;			// Leap seconds
	float UTC_time[6];		// UTC time obtained from the GPS time
	
	/* GPS Position and Velocity */
	// Only the one given by the GPS unit. 
	float gpsPos_ecef[3];	// [gpsX_ecef, gpsY_ecef, gpsZ_ecef]
	float gpsVel_ecef[3];
	float gpsPos_eci[3];	// [gpsX_eci, gpsY_eci, gpsZ_eci]	
	float gpsVel_eci[3];
	
} gps;

extern gps gps_data;

void gps2utc(float GPS_time[2], int leap_sec, float UTC_time[6]);	// TESTED

#endif