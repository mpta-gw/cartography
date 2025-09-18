function [P12, P11, P22, fout] = cohereTimeSeries(d1, d2, fb1, fb2, ...
                                                  dT1, dT2, flow, fhigh);
%
%  cohereTimeSeries --- use to calculate coherence between two time-series
%
%  cohereTimeSeries(d1, d2, fb1, fb2, dT1, dT2, flow, fhigh) calculates
%  the csd and psd's of two (possibly heterodyned) time-series.  Averages 
%  over time intervals, then give the coherence.
%
%  Input:
%
%    d1, d2 = two time-series covering the same length of time
%    fb1, fb2 = two base frequencies associated with the time-series
%    dT1, dT2, = sampling interval of the time series
%    flow, fhigh = frequency span over which to calculate the coherence
%
%  Output:
%
%    P12 = vector with the complex cross power spectral density
%    P11, P22 = two-sided power spectral densities
%    fout = vector of frequencies of the output power spectra
%
%  Routine written by Martin McHugh.
%  Contact mmchugh@loyno.edu
%
%  $Id: cohereTimeSeries.m,v 1.4 2005-04-12 09:09:56 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freso1 = 1/(length(d1)*dT1);
freso2 = 1/(length(d2)*dT2);
if freso1 ~= freso2
    error(['the two data streams must have the same frequency resolution'])
end
if isnan(fb1)  %% check for heterodyning
    f1 = [-1/(2*dT1):freso1:1/(2*dT1)-freso1];
else        
    f1 = [-1/(2*dT1):freso1:1/(2*dT1)-freso1]+fb1;
end

if isnan(fb2)  %% check for heterodyning
    f2 = [-1/(2*dT2):freso2:1/(2*dT2)-freso2];
else        
    f2 = [-1/(2*dT2):freso2:1/(2*dT2)-freso2]+fb2;
end
 
ind1 = find(f1 >= flow & f1 <= fhigh);
ind2 = find(f2 >= flow & f2 <= fhigh);

if floor(ind1/2)~=ind1/2    
    ind1 = ind1(1:length(ind1)-1);  % force the fft to be even number of pts
end
if floor(ind2/2)~=ind2/2    
    ind2 = ind2(1:length(ind2)-1);  % force the fft to be even number of pts
end

w1 = hann(length(d1));
w2 = hann(length(d2));

d1_tild = fftshift(fft(w1.*d1))*dT1;
d2_tild = fftshift(fft(w2.*d2))*dT2;

P12 = ...
sqrt(length(d1)*length(d2))*2*d1_tild(ind1).*conj(d2_tild(ind2))/...
    (sum(abs(w1).^2)*sqrt(length(w2)/length(w1))*length(d1)*dT1);
P11=length(d1)*2*(d1_tild(ind1).*conj(d1_tild(ind1)))/(sum((abs(w1)).^2)*length(d1)*dT1);
P22=length(d2)*2*(d2_tild(ind2).*conj(d2_tild(ind2)))/(sum((abs(w2)).^2)*length(d2)*dT2);
fout = f1(ind1);

return

