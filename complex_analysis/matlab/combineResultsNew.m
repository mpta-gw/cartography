function ...
  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, numSegments] = ...
    combineResultsNew(ccStats, ccSigmas, ccSpectra, sensInts, ...
                      segmentDuration, flow, fhigh, deltaF, ...
                      doOverlap, window1, window2, displayResults)
%
%  
%  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, numSegments] = ...
%    combineResultsNew(ccStats, ccSigmas, ccSpectra, sensInts, ...
%                      segmentDuration, flow, fhigh, deltaF, ...
%                      doOverlap, window1, window2, displayResults)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands.
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
%    numSegments = number of data segments combined
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: combineResultsNew.m,v 1.1 2005-04-07 14:06:07 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract number of segments
numSegments = length(ccStats);

% check consistency of input data 
numFreqs = floor((fhigh-flow)/deltaF)+1;

if numFreqs ~= length(ccSpectra(1,:))
  error('input data has wrong number of frequencies');
end

if numFreqs ~= length(sensInts(1,:))
  error('input data has wrong number of frequencies');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct theoretical variances
ccVars = ccSigmas.^2;

% initialize output data arrays for freq-series
combinedPtEstIntData   = zeros(numFreqs,1);
combinedPtEstIntData_o = zeros(numFreqs,1);
combinedPtEstIntData_e = zeros(numFreqs,1);
combinedSensIntData    = zeros(numFreqs,1);
combinedSensIntData_o  = zeros(numFreqs,1);
combinedSensIntData_e  = zeros(numFreqs,1);

% apply optimal weighting
if ( (doOverlap == false) | (numSegments == 1) )

  ptEstimate = sum((ccStats./ccVars)/sum(1./ccVars)) / segmentDuration;
  errorBar = sqrt(1/sum(1./ccVars)) / segmentDuration;

  for k=1:numFreqs
    combinedPtEstIntData(k,1) = ...
      sum((ccSpectra(:,k)./ccVars)/sum(1./ccVars))/segmentDuration;
  end

  combinedSensIntData = transpose(sum(sensInts,1))*(segmentDuration^2);
 
else

  % construct an array of even and odd indices
  if mod(numSegments,2)==1
    j_o = [1:2:numSegments]';
    j_e = [2:2:numSegments-1]';
  else
    j_o = [1:2:numSegments-1]';
    j_e = [2:2:numSegments]';
  end
  
  % even and odd point estimates
  ptEstimate_o=sum((ccStats(j_o)./ccVars(j_o))/sum(1./ccVars(j_o))) ...
                 / segmentDuration;
  ptEstimate_e=sum((ccStats(j_e)./ccVars(j_e))/sum(1./ccVars(j_e))) ...
                 / segmentDuration;

  % even and odd error bars
  errorBar_o = sqrt( 1/sum(1./ccVars(j_o)) ) / segmentDuration;
  errorBar_e = sqrt( 1/sum(1./ccVars(j_e)) ) / segmentDuration;
  
  % even and odd point estimate integrands 
  for k=1:numFreqs
    combinedPtEstIntData_o(k,1) = ...
      sum( (ccSpectra(j_o,k)./ccVars(j_o) )/sum(1./ccVars(j_o)) ) ...
      / segmentDuration;

    combinedPtEstIntData_e(k,1) = ...
      sum( (ccSpectra(j_e,k)./ccVars(j_e) )/sum(1./ccVars(j_e)) ) ...
      / segmentDuration;
  end

  % even and odd sensitivity integrands 
  combinedSensIntData_o = transpose(sum(sensInts(j_o,:),1))*(segmentDuration^2);
  combinedSensIntData_e = transpose(sum(sensInts(j_e,:),1))*(segmentDuration^2);

  % construct covariance matrix for the even and odd data
  C_oo = errorBar_o^2;
  C_ee = errorBar_e^2;
  
  [w2bar, w4bar, woverlap4bar] = windowFactors(window1, window2);
  sigma2I  = ccVars(1:end-1) / (segmentDuration^2);
  sigma2J  = ccVars(2:end)   / (segmentDuration^2);
  sigma2IJ = 0.5*(woverlap4bar/w4bar)*0.5*(sigma2I+sigma2J);
  
  C_oe = C_oo*C_ee*sum(sigma2IJ./(sigma2I.*sigma2J));
  detC = C_oo*C_ee - C_oe*C_oe;

  % need similar quantities for the sensitivity integrands 
  sensIntsI  = sensInts(1:end-1,:) * (segmentDuration^2);
  sensIntsJ  = sensInts(2:end,:) * (segmentDuration^2);
  sensIntsIJ = 0.5*(woverlap4bar/w4bar)*0.5* ...
                 transpose(sum(sensIntsI,1)+sum(sensIntsJ,1)); 
                                                      
  % construct optimal weighting factors from the covariance matrix
  lambda_o = (C_ee - C_oe)/detC;
  lambda_e = (C_oo - C_oe)/detC;
  
  % final point estimate
  ptEstimate = (ptEstimate_o*lambda_o + ptEstimate_e*lambda_e) ...
               / (lambda_o + lambda_e) ;

  % final error bar 
  errorBar = sqrt(1/(lambda_o + lambda_e));

  % final point estimate integrand 
  combinedPtEstIntData = ...
    (combinedPtEstIntData_o*lambda_o + combinedPtEstIntData_e*lambda_e) ...
    / (lambda_o + lambda_e) ;

  % final sensitivity integrand
  combinedSensIntData = ...
    (combinedSensIntData_o + combinedSensIntData_e - 2*sensIntsIJ) ...
    * ((C_oo*C_ee)/detC);

end

% construct frequency series
combinedPtEstInt = constructFreqSeries(combinedPtEstIntData, flow, deltaF);
combinedSensInt  = constructFreqSeries(combinedSensIntData, flow, deltaF);

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

%fprintf('combineResultsNew.m results ------------------------------------\n');
fprintf('Pt estimate = %e\n', ptEstimate);
fprintf('Error bar   = %e\n', errorBar);
%fprintf('Alternate point estimate calculation = %e\n', ptEstimateAlt);
%fprintf('Alternate error bar calculation = %e\n', errorBarAlt);
fprintf('Number segs = %d\n', numSegments);
%fprintf('combineResultsNew.m results ------------------------------------\n');

return

