function [F]=antenna(det,time,ra,decl,pol)
%  function [F]=antenna(det,time,ra,decl,pol)
%
%  antenna  -- calculates the projection F(ra,decl,pol)
%              on to detector det at sidereal time time
%                    
%  arguments: det  - detector structure containing position r and tensor d
%             time - sidereal time in hours (0..24h)
%             ra   - right ascension in hours of source in the sky
%             decl - right ascension in degrees of source in the sky
%             pol  - polarization angle of GW
%                    relative to north-up coordinate system
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wSidereal=2*pi * (1/365.2425 +1) / 24;
w=2*pi/24;
epz=[1 0 0; 0 -1 0; 0 0 0];
psi=w*(time-ra);
theta=-pi/2+pi/180*decl;
phi=pol;

F=trace(epz*rotateTensor(det.d,phi,theta,psi));

