function [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
          numSegmentsTotal] = ...
  combineResultsFromSingleJob(ccStats, ccSigmas, ccSpectra, sensInts, ...
                              gpsTimes, segmentDuration, ...
                              flow, fhigh, deltaF, ...
                              doOverlap, window1, window2, displayResults)
%
%  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
%   numSegmentsTotal] = ...
%  combineResultsFromSingleJob(ccStats, ccSigmas, ccSpectra, sensInts, ...
%                              gpsTimes, segmentDuration, ...
%                              flow, fhigh, deltaF, ...
%                              doOverlap, window1, window2, displayResults)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands for a single job by calling combineResultsNew 
%  repeatedly.
%
%  Input:
%
%    ccStats = column vector containing the CC statistic values
%    ccSigmas = column vector containing the theoretical sigmas
%    ccSpectra = 2-d array containing complex cc spectra
%    (spectra corresponding to different times are in different rows)
%    sensInts = 2-d array containing real sensitivity integrands
%    (integrands corresponding to different times are in different rows)
%
%    gpsTimes = column vector of GPS start times
%
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    flow, fhigh, deltaF = low, high, and frequency spacing (in Hz)
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%
%    displayResults = 0 don't write results to screen (default = 1)
%
%  Output:
%
%    ptEstimate = optimal point estimate of Omega0
%    errorBar = theoretical error bar for the point estimate
%    combinedPtEstInt = a frequency-series data structure containing the
%      integrand of the optimal point estimate of Omega0
%    combinedSensInt = a frequency-series data structure containing the
%      integrand of 1/ErrorBar^2
%    numSegmentsTotal = total number of data segments combined
%
%  $Id: combineResultsFromSingleJob.m,v 1.3 2005-11-01 04:59:32 charlton Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
ptEstimate = 0;
errorBar = 0;
combinedPtEstInt = constructFreqSeries(0,0,0);
combinedSensInt = constructFreqSeries(0,0,0);
numSegmentsTotal = 0;

% check for non-trivial data
if isempty(gpsTimes) return; end

% combine the data
numSegmentsTotal = length(gpsTimes);
ii=1;
kk=1;

% Eliminate warnings in compiler
ccStats_ovl = [];
ccSigmas_ovl = [];
ccSpectra_ovl = [];
sensInts_ovl = [];

if doOverlap==1

  while ii <= numSegmentsTotal

    clear ccStats_ovl;
    clear ccSigmas_ovl;
    clear ccSpectra_ovl;
    clear sensInts_ovl;

    jj=1;
    ccStats_ovl(jj,1) = ccStats(ii);
    ccSigmas_ovl(jj,1) = ccSigmas(ii);
    ccSpectra_ovl(jj,:) = ccSpectra(ii,:);
    sensInts_ovl(jj,:) = sensInts(ii,:);

    if ii~=numSegmentsTotal
      while gpsTimes(ii+1) == gpsTimes(ii)+segmentDuration/2
        ii=ii+1;
        jj=jj+1;
        ccStats_ovl(jj,1) = ccStats(ii);
        ccSigmas_ovl(jj,1) = ccSigmas(ii);
        ccSpectra_ovl(jj,:) = ccSpectra(ii,:);
        sensInts_ovl(jj,:) = sensInts(ii,:);
        if ii==numSegmentsTotal
          break;
        end  
      end
    end

    [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ignore] = ...
      combineResultsNew(ccStats_ovl, ccSigmas_ovl, ...
                        ccSpectra_ovl, sensInts_ovl, ...
                        segmentDuration, flow, fhigh, deltaF, ...
                        doOverlap, window1, window2, 0);

    % apply appropriate factors of segmentDuration
    ccStats_nonovl(kk,1) = ptEstimate*segmentDuration;
    ccSigmas_nonovl(kk,1) = errorBar*segmentDuration;
    ccSpectra_nonovl(kk,:) = transpose(combinedPtEstInt.data)*segmentDuration;
    sensInts_nonovl(kk,:) = transpose(combinedSensInt.data)/(segmentDuration^2);

    ii=ii+1;
    kk=kk+1;

  end
   
  % combine cc stats, etc. from non-overlapping segments
  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ignore] = ...
    combineResultsNew(ccStats_nonovl, ccSigmas_nonovl, ...
                      ccSpectra_nonovl, sensInts_nonovl, ...
                      segmentDuration, flow, fhigh, deltaF, 0, 0, 0, 0);

else
  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
   numSegmentsTotal] = ...
   combineResultsNew(ccStats, ccSigmas, ccSpectra, sensInts, ...
                     segmentDuration, flow, fhigh, deltaF, 0, 0, 0, 0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write results to the screen
try 
  displayResults;
catch
  displayResults = 1;
end
if displayResults==0
  return
end

ptEstimateAlt = 2*sum(real(combinedPtEstInt.data)*combinedPtEstInt.deltaF);
errorBarAlt = sqrt(1/sum(combinedSensInt.data*combinedSensInt.deltaF));

%fprintf('combineResultsFromSingleJob.m results --------------------------\n');
fprintf('Pt estimate = %e\n', ptEstimate);
fprintf('Error bar   = %e\n', errorBar);
%fprintf('Alternate point estimate calculation = %e\n', ptEstimateAlt);
%fprintf('Alternate error bar calculation = %e\n', errorBarAlt);
fprintf('Number segs = %d\n', numSegmentsTotal);
%fprintf('combineResultsFromSingleJob.m results --------------------------\n');

return
