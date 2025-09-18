function [siderealTime]=GPStoGreenwichMeanSiderealTime(GPS,varargin)
%  function [siderealTime]=GPStoGreenwichMeanSiderealTime(GPS,usemod)
%
%  GPStoGreenwichMeanSiderealTime  -- converts GPS to mean sidereal time
%                                     at greenwich
%                    
%  arguments: GPS          - GPS time in seconds
%             usemod       - apply the "mod 24" if true, default=true
%  output:    siderealTime - mean sidereal time at greenwich in hours
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2
  error('not more than 2 arguments.');
end;
if nargin == 2
  usemod=varargin{1};
else
  usemod=true;
end;

% from International Earth Rotation Service, measured in 2002
% see http://hypertextbook.com/facts/2002/JasonAtkins.shtml
%wearth=72.92115090*1e-6;
wearth=2*pi * (1/365.2425 +1) / 86400;

% same thing in hours / sec
w=wearth/pi*12;  

% GPS time for 1/1/2000 00:00:00 UTC
GPS2000=630720013;
% sidereal time at GPS2000,in hours
% from http://www.csgnetwork.com/siderealjuliantimecalc.html
sT0=6+39/60+51.251406103947375/3600;

siderealTime=w*(GPS-GPS2000)+sT0;
if usemod
  siderealTime=mod(siderealTime,24);
end;
return;
