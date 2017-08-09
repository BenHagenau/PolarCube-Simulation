/*
 * GPS.c
 *
 * Created: 9/27/2013 10:13:18 PM
 *  Author: Alberto
 */ 
#include "GPS.h"
#include "math.h"

/* Function: gps2utc
* Converts the week number and time of week (TOW) to UTC time.
*
* INPUTS:
* GPS_time - GPS Time in the form: 
*					GPS_time[0] = week number (valid numbers: 1-3640 (years 1980-2050))
*					GPS_time[1] = TOW (valid sec values are 0-604799)
* offset   - leap seconds for the GPS time (valid offset values are 0-500)
*
* OUTPUTS:
* UTC_time - matrix in the form
*				UTC_time[0] = year (yyyy) 
*				UTC_time[1] = month
*				UTC_time[2] = day
*				UTC_time[3] = hour (24hr)
*				UTC_time[4] = min
*				UTC_time[5] = sec
*/
void gps2utc(float GPS_time[2], int leap_sec, float UTC_time[6])
{
	float GPS_week;
	float GPS_sec;
	float UTC_temp[6];
	float gpsday;
	float gpssec;
	float total_days;
	int   leap_year, I_leap, I_no_leap;
	float day_of_year;

	GPS_week = GPS_time[0]; 
	GPS_sec  = GPS_time[1];

	// Compute gpsday and gps seconds since start of GPS time  
	gpsday = GPS_week * 7 + GPS_sec / 86400; 
	gpssec = GPS_week * 7 * 86400 + GPS_sec; 

	// get the integer number of days 
	total_days = (float) floor(gpsday); 

	/* temp is the number of completed years since the last leap year (0-3) 
	% the calculation begins by computing the number of full days since 
	% the beginning of the last leap year.  This is accomplished through 
	% the rem statement.  Since GPS time started at 
	% 00:00 on 6 January 1980, five days must be added to the total number 
	% of days to ensure that the calculation begins at the beginning of a 
	% leap year.  By subtracting one from this result, the extra day in 
	% the first year is effectively removed, and the calculation can 
	% simply be computed by determining the number of times 365 divides 
	% into the number of days since the last leap year.  On the first day 
	% of a leap year, the result of this calculation is -1 
	% so the second statement is used to trap this case. */
 
	float temp = (float)floor(((float)fmod((total_days+5),1461)-1) / 365); 

	temp = (float)floor(((float)fmod((total_days+5),1461)-1) / 365); 
	if (temp < 0)
	{
	  temp = 0;  
	}
 
	// compute the year 
	UTC_time[0] = 1980 + 4 * (float)floor((total_days + 5) / 1461) + temp;

	/* data matrix with the number of days per month for searching  
	% for the month and day 
	% days in full months for leap year */
	float leapdays[] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};   
	// days in full months for standard year 
	float noleapdays[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365}; 

	// Leap year flag 
	// determine which input years are leap years 
	leap_year = (float)fmod((UTC_time[0]-1980),4);
	if (leap_year == 0)
	{
		I_leap    = 1;     // leap year
		I_no_leap = 0;
	}
	else
	{
		I_leap    = 0;
		I_no_leap = 1;  // standard year
	}

	// establish the number of days into the current year 
	// leap year 
	if (I_leap)
	{
	  day_of_year = (float)fmod((total_days + 5),1461) + 1;                         
	}
 
	// standard year 
	if (I_no_leap) 
	{
	  day_of_year = (float)fmod((float)fmod((total_days + 5),1461) - 366, 365) + 1;   
	}

	// generate the month, loop over the months 1-12 and separate out leap years 
	for (int iii = 0; iii < 13; iii++)
	{
	  if (I_leap )
	  {
		if (day_of_year > leapdays[iii])
		{
			UTC_time[1] = iii + 1; 
		}
	  }
   
	  if (I_no_leap)
	  {
		if (day_of_year > noleapdays[iii])
		{
			UTC_time[1] = iii + 1;
		}
	  }
	}

	// use the month and the matrix with days per month to compute the day  
	if (I_leap)
	{
		int month = UTC_time[1] - 1;
		UTC_time[2] = day_of_year - leapdays[9]; 
	}
 
	if (I_no_leap)
	{
		int month = UTC_time[1] - 1;
		UTC_time[2] = day_of_year - noleapdays[month]; 
	}
 
	// compute the hours 
	float fracday = (float)fmod(GPS_sec, 86400);              // in seconds! 
 
	UTC_time[3] = (float)floor(fracday / 86400 * 24); 
 
	// compute the minutes  
	UTC_time[4] = (float)floor((fracday - UTC_time[3] * 3600) / 60 ); 
 
	// compute the seconds 
	UTC_time[5] = fracday - UTC_time[3] * 3600 - UTC_time[4] * 60; 

	/* Compensate for leap seconds */ 
	UTC_time[5] = UTC_time[5] - leap_sec; 
 
	// Check to see if leap_sec offset causes a negative number of seconds 
	if (UTC_time[5] < 0)
	{
		UTC_time[4] = UTC_time[4] - 1; 
		UTC_time[5] = UTC_time[5] + 60; 
	}

	// Check to see if the leap second offset causes a negative number of minutes 
	if (UTC_time[4] < 0)
	{
		UTC_time[3] = UTC_time[3] - 1; 
		UTC_time[4] = UTC_time[4] + 60; 
	}
 
	// Check to see if the leap second offset causes a negative number of hours 
	if (UTC_time[3] < 0)
	{
		UTC_time[2] = UTC_time[2] - 1; 
		UTC_time[3] = UTC_time[3] + 24; 
	}

	// Check to see if this causes a 0 day value 
	if (UTC_time[2] <= 0)
	{
		UTC_time[1] = UTC_time[1] - 1; 
	}
	if (UTC_time[1] <= 0)
	{
		UTC_time[0] = UTC_time[0] - 1; 
		UTC_time[1] = UTC_time[1] + 12; 
	}
 
	/* Leap year flag 
	 % determine which input years are leap years */
	leap_year = (float)fmod((UTC_time[0]-1980),4);
	if (leap_year == 0)
	{
		I_leap    = 1;     // leap years
		I_no_leap = 0;
	}
	else
	{
		I_leap    = 0;
		I_no_leap = 1;  // standard years
	}
 
	if (I_leap)
	{
		int month   = UTC_time[1] - 1;
		UTC_time[2] = leapdays[month + 1] - leapdays[month]; 
	}

	if (I_no_leap)
	{
		if (UTC_time[2] <= 0)
		{
			int month   = UTC_time[1] - 1;
			UTC_time[2] = noleapdays[month + 1] - noleapdays[month]; 
		}
	}
}