function y = sumTerms(x, jlow, jhigh)
%
%  sumTerms --- a routine needed for coarse-graining
%
%  sumTerms(x, jlow, jhigh) returns an array of values corresponding to
%  the array x[j] summed from jlow to jhigh:
%
%  y[k] = sum_{j=jlow[k]}^{jhigh[k]} x[j] 
%
%  where jlow[k] <= jhigh[k] (k = 1, 2, ... Ny < Nx)
%
%  NOTE: This routine should eventually be turned into a MEX function
%  to improve the efficiency of coarseGrain.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: sumTerms.m,v 1.2 2005-02-24 15:11:32 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ny = length(jlow);
y = zeros(Ny,1);

for k=1:1:Ny
  y(k) = sum( x(jlow(k):jhigh(k)) );
end

return

