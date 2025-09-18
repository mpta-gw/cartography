function minmax = varLimit(valueIn, highLow)
%
%  varLimit(valueIn, highLow) takes a continuous input valueIn and 
%  determines the closest integer = +/-(1 2 3 5 8 10)*10^n for a plot 
%  with proper accounting of whether valueIn is > or < 0. 
%  Use highLow = 1 for setting an upper axis limit and 
%      highlow = 0 for a lower axis limit; 
%  NOTE all  values ~= 0 are set to 1.
%
%  Routine written by Albert Lazzarini.
%  Contact lazz@ligo.caltech.edu
%
%  $Id: varLimit.m,v 1.1 2005-10-19 17:51:01 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

direction = sign(valueIn);
value = abs(valueIn);
uplow = highLow;

if uplow ~= 0
  uplow = 1;
end

if direction < 0
  uplow = 1 - uplow;
end

% calculate quantized limits for plots or histograms
nearestTen = log10([1 2 3 5 8 10]);
lnym = log10(value); 
flnym = floor(lnym); 

if uplow == 1
  inds = find(lnym<flnym+nearestTen);
  minmax = direction*10^(flnym+nearestTen(inds(1)));
else
  inds = find(lnym>flnym+nearestTen);
  minmax = direction*10^(flnym+nearestTen(inds(end)));
end

return

