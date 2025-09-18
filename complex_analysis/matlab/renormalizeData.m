function [ccStatNew, ccSigmaNew, ccSpectrumNew, sensIntNew] = ...
  renormalizeData(ccSpectrum, sensInt, modifyFilter, flow, fhigh, deltaF)
%
%  renormalizeData -- renormalize CC statistic value, spectrum, etc. 
%  for a different choice of optimal filter
%
%  [ccStatNew, ccSigmaNew, ccSpectrumNew, sensIntNew] = ...
%  renormalizeData(ccSpectrum, sensInt, modifyFilter, flow, fhigh, deltaF)
%  returns a renormalized CC statistic value, theoretical sigma, CC 
%  spectrum, and sensitivity integrand for a different choice of optimal
%  filter.  Useful for reanalysing data when applying additional frequency 
%  masking or for different spectral indices, etc. 
%
%  Input:
%
%    ccSpectrum = row vector containing the values of the CC spectrum
%    sensInt = row vector containing the values of the sensitivity 
%      integrand (i.e., the integrand of 1/sigma^2)
%    modifyFilter = row vector containing values for modifying the 
%      optimal filter function (e.g., if modifyFilter = all 1's, then 
%      there will be no change to optimal filter)
%    flow, fhigh, deltaF = low, high, and frequency spacing (in Hz)
%    
%  Output:
%
%    ccStatNew = renormalized CC statistic value and
%    ccSigmaNew = renormalized theoretical sigma
%    ccSpectrumNew = row vector containing the values of the renormalized
%      CC spectrum
%    sensIntNew = row vector containing the values of the renormalized
%      sensitivity integrand
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: renormalizeData.m,v 1.1 2005-04-07 14:01:31 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test for valid freq-series data
numFreqs = floor((fhigh-flow)/deltaF)+1;

if numFreqs ~= length(ccSpectrum)
  error('cc spectrum has wrong number of frequencies');
end
if numFreqs ~= length(sensInt)
  error('sensitivity integrand has wrong number of frequencies');
end
if numFreqs ~= length(modifyFilter)
  error('modify filter function has wrong number of frequencies');
end
                                                                               
% modify sensitivity integrand for different filter and calculate new 
% cc sigma 
sensIntNew = sensInt.*(modifyFilter.^2);
ccSigmaNew = sqrt(1/sum(sensIntNew*deltaF));
                                                                                
% modify cc spectrum for different filter and calculate new cc stat
ccSigmaOld = sqrt(1/sum(sensInt*deltaF));
factor = (ccSigmaNew/ccSigmaOld)^2;
ccSpectrumNew = factor * (ccSpectrum.*modifyFilter);
ccStatNew = 2*deltaF*sum(real(ccSpectrumNew));
   
return
