function [ptEstIntegrandNew, sensIntNew] = ...
  diffFreqMask(ptEstIntegrand, sensInt, freqsToRemove, nBinsToRemove, ...
               outputFileNamePrefix, displayResults, ...
               figureNumber, figureLegend, complexFlag);
%
%  [ptEstIntegrandNew, sensIntNew] = ...
%  diffFreqMask(ptEstIntegrand, sensInt, freqsToRemove, nBinsToRemove, ...
%               outputFileNamePrefix, displayResults, ...
%               figureNumber, figureLegend);
% 
%  calculate point estimate and error bar when different frequency masking 
%  is used
%
%  Input:
% 
%    ptEstIntegrand = a freq-series data structure containing the
%      integrand of the point estimate of Omega
%    sensInt = a freq-series data structure containing the sensitivity
%      integrand (i.e., the integrand of 1/errorbar^2)
%    freqsToRemove = column vector containing freq values (in Hz) that
%      will be set to zero
%    nBinsToRemove = column vector containing the number of bins to 
%      remove around each of the discrete freqs selected to be removed
%    outputFileNamePrefix = prefix of filenames that will contain
%      the new point estimate and sensitivity integrands;
%      prefix of filenames for figures if displaying results.
%    displayResults = 0 don't display results to screen (default = 1)
%    figureNumber = number of first figure (if displayResults = 1)
%    figureLegend = text for figure legend (if displayResults = 1)
%
%  Output:
%
%    ptEstIntegrandNew: a freq-series data structure containing the
%      integrand of the point estimate of Omega appropriately modified 
%      for the different frequency masking 
%    sensIntNew: a freq-series data structure containing the sensitivity 
%      integrand appropriately modified for the different frequency masking
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%  
%  $Id: diffFreqMask.m,v 1.15 2005-11-16 20:37:21 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  complexFlag;
catch
  complexFlag = false;
end

ddmmyyyyhhmmss  = datestr(now);

% extract frequency series information
deltaF = ptEstIntegrand.deltaF;
numFreqs = length(ptEstIntegrand.data);
flow = ptEstIntegrand.flow;
freqs = flow+deltaF*[0:numFreqs-1]';
fhigh = freqs(end);                                                                                
% calculate original point estimate and error bar
if complexFlag == false
  pointEstimate = 2*deltaF*sum(real(ptEstIntegrand.data));
else
  pointEstimate = deltaF*sum(ptEstIntegrand.data);
end;  
errorBar = sqrt(1/sum(sensInt.data*deltaF));
                                                       
% construct freq mask
data = constructFreqMask(flow, fhigh, deltaF, ...
                         freqsToRemove, nBinsToRemove, 1);
mask = constructFreqSeries(data, flow, deltaF);

% modify sens int for different freq mask and calculate new error bar
data = sensInt.data .* mask.data;
sensIntNew = constructFreqSeries(data, flow, deltaF);
errorBarNew = sqrt(1/sum(sensIntNew.data*deltaF));

% modify cc spectrum for different freq mask and calculate new pt estimate
factor = (errorBarNew/errorBar)^2; 
data = factor * (ptEstIntegrand.data .* mask.data);
ptEstIntegrandNew = constructFreqSeries(data, flow, deltaF);
pointEstimateNew = 2*deltaF*sum(real(ptEstIntegrandNew.data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write modified cc spectrum, sens integrand to file
                                                                                
% cc spectrum -----------------------------------------------------------
filename = [outputFileNamePrefix '_diffFreqMask_ptEstIntegrand.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\t%s\t%s\n', ...
        '%error bar', 'freq (Hz)', 'pt est int (re)', 'pt est int (im)');
for ii=1:numFreqs
  fprintf(fid, '%e\t%e\t%e\t%e\n', errorBarNew, freqs(ii), ...
          real(ptEstIntegrandNew.data(ii)), imag(ptEstIntegrandNew.data(ii)));
end
fclose(fid);
                                                                                
% sens integrand --------------------------------------------------------
filename = [outputFileNamePrefix '_diffFreqMask_sensIntegrand.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\t%s\n', ...
        '%error bar', 'freq (Hz)', 'sens integrand');
for ii=1:numFreqs
  fprintf(fid, '%e\t%e\t%e\n', errorBarNew, freqs(ii), sensIntNew.data(ii));
end
fclose(fid);
                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display final results if desired
try
  displayResults;
catch
  displayResults = 1;
end
if displayResults==0
  return
end
                                                                                
% display original results
fprintf('ORIGINAL RESULTS:\n');
fprintf('Pt estimate = %e\n', pointEstimate);
fprintf('Error bar   = %e\n', errorBar);

% display renormalized results
fprintf('RESULTS FOR DIFFERENT FREQ MASKING:\n');
fprintf('Pt Estimate = %e\n', pointEstimateNew);
fprintf('Error bar   = %e\n', errorBarNew);

% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                
% modified cc spectrum
figure(figureNumber)
plot(freqs, real(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Real(Integrand of point estimate)', 'FontSize', 14);
title('Modified point estimate integrand (diff freq masking)','fontsize',14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffFreqMask_ptEstIntegrand_real'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+1)
plot(freqs, imag(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Imag(Integrand of point estimate)', 'FontSize', 14);
title('Modified point estimate integrand (diff freq masking)','fontsize',14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffFreqMask_ptEstIntegrand_imag'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+2)
plot(freqs, abs(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Abs(Integrand of point estimate)', 'FontSize', 14);
title('Modified point estimate integrand (diff freq masking)','fontsize',14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffFreqMask_ptEstIntegrand_abs'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+3)
cumpointest = cumsum(2*real(ptEstIntegrandNew.data)*deltaF);
plot(freqs, cumpointest);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative point estimate', 'FontSize', 14);
title('Modified cumulative point estimate (diff freq masking)','fontsize',14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffFreqMask_ptEstIntegrand_cum'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
% modified sensitivity integrand
figure(figureNumber+4)
plot(freqs, sensIntNew.data);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Sensitivity Integrand', 'FontSize', 14);
title('Modified sensitivity integrand (diff freq masking)','fontsize',14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffFreqMask_sensIntegrand'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+5)
cumsens = cumsum(sensIntNew.data);
cumsens = cumsens/cumsens(end);
plot(freqs, cumsens);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative sensitivity', 'FontSize', 14);
title('Modified cumulative sensitivity (diff freq masking)','fontsize',14);
legend(figureLegend, 2);
filename = [outputFileNamePrefix '_diffFreqMask_sensIntegrand_cum'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');

return
