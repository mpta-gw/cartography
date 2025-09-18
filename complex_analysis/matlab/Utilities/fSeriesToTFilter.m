function timeSeries = fSeriesToTFilter(freqSeries, N, deltaT, inferNegFreqs)
%
%  fSeriesToTFilter -- Generate the time-domain filter corresponding to
%                      a frequency series
%
%  fSeriesToTFilter(freqSeries, N, deltaT) generates from the frequency
%  series freqSeries the corresponding time-domain filter with N points
%  spaced deltaT apart
%
%  This also includes undoing a previous coarse-graining step to get the
%  frequency-domain resolution consistent
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: fSeriesToTFilter.m,v 1.3 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derive medatata for fine-grained frequency series
deltaF = 1/(N*deltaT);
flow = -floor(N/2)*deltaF;

fineGrained = reverseCoarseGrain(freqSeries, flow, deltaF, N, inferNegFreqs);

timeSeries.deltaT = deltaT;
timeSeries.tlow = -floor(N/2)*deltaT;
timeSeries.data = ifftshift(ifft(fftshift(fineGrained.data)))/deltaT;

if inferNegFreqs
  timeSeries.data = real(timeSeries.data);
end

return
