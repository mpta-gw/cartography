function [out]=rotateVector(in,phi,theta,psi)
%  function [out]=rotateVector(in,phi,theta,psi)
%
%  rotateVector   -- rotate a vector by
%                     psi around z-axis, then by
%                     theta around y-axis, then by
%                     phi around z-axis
%                    
%
%  rotateVector(in,phi,theta,psi) rotates the 3d vector in by
%  psi around z-axis, then by theta around y-axis, then by
%  phi around z-axis
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                             
r1=[ cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
r2=[ cos(theta) 0 sin(theta) ; 0 1 0; -sin(theta) 0 cos(theta)];
r3=[ cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];

r=r1*r2*r3;

out=r*in;

