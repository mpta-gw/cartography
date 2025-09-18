function [out]=rotateTensor(in,phi,theta,psi)
%  function [out]=rotateTensor(in,phi,theta,psi)
%
%  rotateTensor   -- rotate a tensor by
%                     psi around z-axis, then by
%                     theta around y-axis, then by
%                     phi around z-axis
%                    
%
%  rotateTensor(in,phi,theta,psi) rotates the 3x3 tensor in by
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

out=r*in*r';

