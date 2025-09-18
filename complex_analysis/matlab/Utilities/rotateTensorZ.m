function [out]=rotateTensorZ(in,phi)
%  function [out]=rotateTensorZ(in,phi)
%
%  rotateTensorZ   -- rotate a tensor by
%                     phi around z-axis
%                    
%
%  rotateTensorZ(in,phi) rotates the 3x3 tensor in by
%  phi around z-axis
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=[ cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

out=r*in*r';

