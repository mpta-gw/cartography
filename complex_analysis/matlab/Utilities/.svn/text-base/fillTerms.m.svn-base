function x = fillTerms(y, jlow, jhigh, Nx)
%
%  fillTerms --- a routine needed for reverse coarse-graining
%
%  fillTerms(y, jlow, jhigh) returns an array of values copied from
%  the array y for jlow to jhigh:
%
%  x(jlow[1]:jhigh[1]) = y[1]
%  x(jlow[2]:jhigh[2]) = y[2]
%    ...
%  x(jlow[Ny]:jhigh[Ny]) = y[Ny]
%
%  where jlow[k] <= jhigh[k] (k = 1, 2, ... Ny < Nx)
%
%  NOTE: This routine should eventually be turned into a MEX function
%  to improve the efficiency of reverseCoarseGrain.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: fillTerms.m,v 1.3 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = zeros(Nx,1);

Ny = length(y);

for k=1:1:Ny
  x(jlow(k):jhigh(k)) = y(k);
end

return

