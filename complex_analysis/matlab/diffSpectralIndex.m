function [ptEstIntegrandNew, sensIntNew] = ...
  diffSpectralIndex(ptEstIntegrand, sensInt, alpha, fRef, ...
                    outputFileNamePrefix, displayresults, ...
                    figureNumber, figureLegend, complexFlag)
%
%  [ptEstIntegrandNew, sensIntNew] = ...
%  diffSpectralIndex(ptEstIntegrand, sensInt, alpha, fRef, ...
%                    outputFileNamePrefix, displayresults, ...
%                    figureNumber, figureLegend, complexFlag)
% 
%  calculate point estimate and error bar for different spectral index 
%  alpha
%
%  Input:
%
%    ptEstIntegrand = a freq-series data structure containing the
%      integrand of the point estimate of Omega
%    sensInt = a freq-series data structure containing the sensitivity
%      integrand (i.e., the integrand of 1/errorbar^2)
%    alpha = spectral index (Omega(f) = OmegaRef*(f/fRef)^alpha)
%    fRef = reference frequency
%    outputFileNamePrefix = prefix of filenames that will contain
%      the new point estimate and sensitivity integrands;
%      prefix of filenames for figures if displaying results.
%    displayresults = 0 don't display results to screen (default = 1)
%    figureNumber = number of first figure (if displayresults = 1)
%    figureLegend = text for figure legend (if displayresults = 1)
%
%  Output:
%
%    ptEstIntegrandNew: a freq-series data structure containing the
%      integrand of the point estimate of Omega appropriately modified 
%      for the different spectral index
%    sensIntNew: a freq-series data structure containing the sensitivity 
%      integrand appropriately modified for the different spectral index
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: diffSpectralIndex.m,v 1.19 2006-04-14 23:39:22 whelan Exp $
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

% calculate original point estimate and error bar
if complexFlag == false
  pointEstimate = 2*deltaF*sum(real(ptEstIntegrand.data));
else
  pointEstimate = deltaF*sum(ptEstIntegrand.data);
end;  
errorBar = sqrt(1/sum(sensInt.data*deltaF));

% modify sens int for different power law and calculate new error bar
data = sensInt.data .* (freqs/fRef).^(2*alpha);
sensIntNew = constructFreqSeries(data, flow, deltaF);
errorBarNew = sqrt(1/sum(sensIntNew.data*deltaF));

% modify cc spectrum for different power law and calculate new pt estimate
factor = (errorBarNew/errorBar)^2; 
data = factor * (ptEstIntegrand.data .* (freqs/fRef).^alpha);
ptEstIntegrandNew = constructFreqSeries(data, flow, deltaF);
pointEstimateNew = 2*deltaF*sum(real(ptEstIntegrandNew.data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write modified cc spectrum, sens integrand to file
 
% cc spectrum -----------------------------------------------------------
filename = [outputFileNamePrefix '_diffSpectralIndex_ptEstIntegrand.dat'];
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
filename = [outputFileNamePrefix '_diffSpectralIndex_sensIntegrand.dat'];
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
  displayresults;
catch
  displayresults = 1;
end
if displayresults==0
  return
end
                                                                                
% display original results
fprintf('ORIGINAL RESULTS:\n');
fprintf('Pt estimate = %e\n', pointEstimate);
fprintf('Error bar   = %e\n', errorBar);

% display renormalized results
fprintf('RESULTS FOR DIFFERENT SPECTRAL INDEX ALPHA:\n');
fprintf('Spectral index  = %d\n', alpha);
fprintf('Pt Estimate = %e\n', pointEstimateNew);
fprintf('Error bar   = %e\n', errorBarNew);

% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modified cc spectrum
figure(figureNumber)
plot(freqs, real(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Real(Integrand of point estimate)', 'FontSize', 14);
titletext = ['Modified point estimate integrand (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffSpectralIndex_ptEstIntegrand_real'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+1)
plot(freqs, imag(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Imag(Integrand of point estimate)', 'FontSize', 14);
titletext = ['Modified point estimate integrand (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffSpectralIndex_ptEstIntegrand_imag'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                               
figure(figureNumber+2)
plot(freqs, abs(ptEstIntegrandNew.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Abs(Integrand of point estimate)', 'FontSize', 14);
titletext = ['Modified point estimate integrand (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffSpectralIndex_ptEstIntegrand_abs'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+3)
cumpointest = cumsum(2*real(ptEstIntegrandNew.data)*deltaF);
plot(freqs, cumpointest);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative point estimate', 'FontSize', 14);
titletext = ['Modified cumulative point estimate (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffSpectralIndex_ptEstIntegrand_cum'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
% modified sensitivity integrand 
figure(figureNumber+4)
plot(freqs, sensIntNew.data);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Sensitivity Integrand', 'FontSize', 14);
titletext = ['Modified sensitivity integrand (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend);
filename = [outputFileNamePrefix '_diffSpectralIndex_sensIntegrand'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
figure(figureNumber+5)
cumsens = cumsum(sensIntNew.data);
cumsens = cumsens/cumsens(end);
plot(freqs, cumsens);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative sensitivity', 'FontSize', 14);
titletext = ['Modified cumulative sensitivity (alpha=' num2str(alpha) ', fRef=' num2str(fRef) ' Hz)'];
title(titletext, 'fontsize', 14);
legend(figureLegend, 2);
filename = [outputFileNamePrefix '_diffSpectralIndex_sensIntegrand_cum'];
print(gcf, '-depsc2', filename);
saveas(gcf, filename, 'fig');
                                                                                
return
