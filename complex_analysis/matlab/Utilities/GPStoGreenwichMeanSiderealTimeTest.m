function [siderealTime]=GPStoGreenwichMeanSiderealTimeTest(GPS)
%  function [siderealTime]=GPStoGreenwichMeanSiderealTime(GPS)
%
%  GPStoGreenwichMeanSiderealTimeTest  -- converts GPS to mean sidereal time
%                                     at greenwich
%                    
%  arguments: GPS - GPS time in seconds
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   /* cf. S. Aoki et al., A&A 105, 359 (1982) eqs. 13 & 19 */
%   /* also cf. http://aa.usno.navy.mil */
%   /* Note: 00h UT 01 Jan 2000 has JD=2451544.5 and GPS=630720013 */
   JD_12h_01_Jan_2000     = 2451545.0;
   JD_00h_01_Jan_2000     = 2451544.5;
   GPS_00h_01_Jan_2000    = 630720013;
   TAIUTC_00h_01_Jan_2000 = 32;% /* leap seconds: TAI - UTC */

%   REAL8 t;
%   REAL8 dpU;
%   REAL8 TpU;
%   REAL8 gmst;
gpssec=GPS;
gpsnan=0;
taiutc=32;

%   /* compute number of seconds since 00h UT 01 Jan 2000 */
   t  = gpssec - GPS_00h_01_Jan_2000;
   t =t+ 1e-9 * gpsnan;
   t =t+ taiutc - TAIUTC_00h_01_Jan_2000;

%   /* compute number of days since 12h UT 01 Jan 2000 */
   dpU  = floor( t / ( 24.0 * 3600.0 ) ); %/* full days since 0h UT 01 Jan 2000 */
   dpU =dpU+ JD_00h_01_Jan_2000 - JD_12h_01_Jan_2000;% /* i.e., -0.5 */

%   /* compute number of centuries since 12h UT 31 Dec 1899 */
   TpU = dpU / 36525.0;

%   /* compute the gmst at 0h of the current day */
   gmst = 24110.54841 ...
     + TpU * ( 8640184.812866 ...
         + TpU * ( 0.093104 ...
           - TpU * 6.2e-6 ) ); %/* seconds */

%   /* add the sidereal time since the start of the day */
   t = mod( t, 24.0 * 3600.0 ); %/* seconds since start of day */
   gmst =gmst+ t * 1.002737909350795; %/* corrections omitted */

%   /* convert to fractions of a day and to radians */
   gmst = mod( gmst / ( 24.0 * 3600.0 ), 1.0 ); %/* fraction of day */
%   /* gmst *= 2.0 * M_PI; */ /* radians */ /* 10/15/04 gam */
%   gmst =gmst* 2*pi; /* radians */
   
%   return gmst;
siderealTime = gmst *24; % in hours