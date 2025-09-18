function [x, index1, index2, frac1, frac2] ...
      = reverseCoarseGrain(y, flowx, deltaFx, Nx, inferNegFreqs)
%
%  reverseCoarseGrain --- fine-grain a frequency-series
%
%  reverseCoarseGrain(y, flowx, deltaFx, Nx, inferNegFreqs) returns a 
%  fine-grained frequency-series structure from the coarse-grained 
%  frequency-series y.  inferNegFreqs is a flag which determines whether 
%  negative frequencies are also set to make the Fourier transform of a 
%  real time series.
%
%  reverseCoarseGrain also returns the indices of the lowest and highest 
%  frequency bins of x that overlap with y (0 <= index1 <= length(x);
%  1 <= index2 <= length(x)+1) and the fractional contributions from these
%  frequency bins.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: reverseCoarseGrain.m,v 1.3 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract metadata
Ny = length(y.data);
flowy = y.flow;
deltaFy = y.deltaF;

% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that frequency spacings and number of elements are > 0
if ( (deltaFx <= 0) | (deltaFy <= 0) | (Ny <= 0) )
  error('bad input arguments\n');
end;

% check that freq resolution of x is finer than desired freq resolution
if ( deltaFy < deltaFx )
  error('deltaF coarse-grain < deltaF fine-grain\n');
end;

% check desired start frequency for coarse-grained series
if ( (flowy - 0.5*deltaFy) < (flowx - 0.5*deltaFx) )
  error('desired coarse-grained start frequency is too low\n');
end;

% check desired stop frequency for coarse-grained series
fhighx = flowx + (Nx-1)*deltaFx;
fhighy = flowy + (Ny-1)*deltaFy;

if ( (fhighy + 0.5*deltaFy) > (fhighx + 0.5*deltaFx) )
  error('desired coarse-grained stop frequency is too high\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indices for coarse-grained series
i = transpose(0:1:Ny-1);

% calculate the low and high indices of fine-grained series 
jlow  = 1 + floor( (flowy + (i-0.5)*deltaFy - flowx - 0.5*deltaFx)/deltaFx ); 
jhigh = 1 + floor( (flowy + (i+0.5)*deltaFy - flowx - 0.5*deltaFx)/deltaFx ); 

if (jhigh(Ny) > Nx-1)
  error('not enough data in coarse-grained series')
end;

index1 = jlow(1);
index2 = jhigh(end);

% calculate fractional contributions
fraclow  = (flowx + (jlow+0.5)*deltaFx - flowy - (i-0.5)*deltaFy)/deltaFx;
frachigh = (flowy + (i+0.5)*deltaFy - flowx - (jhigh-0.5)*deltaFx)/deltaFx;
  
frac1 = fraclow(1);
frac2 = frachigh(end);

% fill in fine-grained values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x.data = zeros(Nx,1);
x.flow = flowx;
x.deltaF = deltaFx;

% fill in middle terms
% NOTE:  since the corresponding operation in coarseGrain.m takes more
% time than anything else, it's in a separate routine called \
% fillTerms() which should probably be converted into a MEX file
jtemp = jlow+2;
x.data = fillTerms(y.data, jtemp, jhigh, Nx);

% fill in edges
x.data(jlow(1:Ny)+1) = fraclow(1:Ny) .* y.data;
x.data(jhigh(1:Ny)+1) = x.data(jhigh(1:Ny)+1) + frachigh(1:Ny) .* y.data;

% if desired, include reverse-coarse-graining of implied negative frequencies
if inferNegFreqs
  yflipped.flow = - ( flowy + (Ny-1)*deltaFy );
  yflipped.deltaF = y.deltaF;
  yflipped.data = flipud(y.data);
  xflipped = reverseCoarseGrain(yflipped, flowx, deltaFx, Nx, false);
  x.data = x.data + xflipped.data;
end

return
