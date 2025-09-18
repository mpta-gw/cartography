function [combinedData, errorBar, numSegments] = ...
  combineResults(combineWhat, inputData, ccSigmas, segmentDuration, ...
                 flow, fhigh, deltaF, doOverlap, window1, window2)
% 
% [combinedData, errorBar, numSegments] = ...
%   combineResults(combineWhat, inputData, ccSigmas, segmentDuration, ...
%                  flow, fhigh, deltaF, doOverlap, window1, window2)
%
%  Optimally combine cc stats, cc spectra, or sensitivity integrands
%  depending on input arguments 
%
%  Input:
%
%    combineWhat = 'stats', 'spectra', 'sensints'
%
%    inputData: 
%     -for cc stats: a column vector containing the cc statistic values
%     -for cc spectra: 2-d array containing complex cc spectra 
%      (spectra corresponding to different times are in different rows)
%     -for sens ints: 2-d array containing real sensitivity integrands
%      (integrands corresponding to different times are in different rows)
%
%    ccSigmas = column vector of theoretical sigmas
%
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    flow, fhigh, deltaF: min, max, and frequency spacing needed for 
%    combining cc spectra and sensitivity integrands 
%    (note: we check that floor((fhigh-flow)/deltaF)-1 is consistent 
%    with the number of frequencies given by the columns of inputData)
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series 
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%
%  Output:
%
%    combinedData:
%     -for cc stats: the optimal point estimate of Omega0
%     -for cc spectra: a frequency-series data structure containing the 
%      integrand of the optimal point estimate of Omega0
%     -for sens ints: a frequency-series data structure containing the 
%      integrand of 1/ErrorBar^2
%
%    errorBar = theoretical error bar for the point estimate
%
%    numSegments = number of data segments combined
%  
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: combineResults.m,v 1.12 2005-02-25 09:36:30 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract information
numSegments = length(ccSigmas);

% check consistency of input data 
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) | ...
     strncmp(combineWhat, 'sensints', length(combineWhat)) )

  % test valid for freq-series data
  numFreqs = floor((fhigh-flow)/deltaF)+1;
  if numFreqs ~= length(inputData(1,:))
    error('input data has wrong number of frequencies');
  end

else
  % test valid for any input data
  if numSegments ~= length(inputData(:,1))
    error('input data has wrong number of segments');
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output data arrays for freq-series
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) | ...
     strncmp(combineWhat, 'sensints', length(combineWhat)) )
  data   = zeros(numFreqs,1);
  data_o = zeros(numFreqs,1);
  data_e = zeros(numFreqs,1);
end

% construct theoretical variances
ccVars = ccSigmas.^2;

% apply optimal weighting
if ( (doOverlap == false) | (numSegments == 1) )

  % error bar -------------------------------------------------------------
  errorBar = sqrt(1/sum(1./ccVars)) / segmentDuration;

  % cc stats --------------------------------------------------------------
  if ( strncmp(combineWhat, 'stats', length(combineWhat)) )
    combinedData = sum((inputData./ccVars)/sum(1./ccVars)) ...
                   / segmentDuration;
  end

  % cc spectra ------------------------------------------------------------
  if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
    for k=1:numFreqs
      data(k,1) = sum((inputData(:,k)./ccVars)/sum(1./ccVars)) ...
                  / segmentDuration;
    end
  end

  % sens ints ------------------------------------------------------------
  if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
    data = transpose(sum(inputData,1)) * (segmentDuration^2);
  end
 
else

  % construct an array of even and odd indices
  if mod(numSegments,2)==1
    j_o = [1:2:numSegments]';
    j_e = [2:2:numSegments-1]';
  else
    j_o = [1:2:numSegments-1]';
    j_e = [2:2:numSegments]';
  end
  
  % error bars ------------------------------------------------------------
  errorBar_o = sqrt( 1/sum(1./ccVars(j_o)) ) / segmentDuration;
  errorBar_e = sqrt( 1/sum(1./ccVars(j_e)) ) / segmentDuration;
  
  % cc stats --------------------------------------------------------------
  if ( strncmp(combineWhat, 'stats', length(combineWhat)) )
    combinedData_o=sum((inputData(j_o)./ccVars(j_o))/sum(1./ccVars(j_o))) ...
                   / segmentDuration;
    combinedData_e=sum((inputData(j_e)./ccVars(j_e))/sum(1./ccVars(j_e))) ...
                   / segmentDuration;
  end 

  % cc spectra ------------------------------------------------------------
  if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
    for k=1:numFreqs
      data_o(k,1) = ...
        sum( (inputData(j_o,k)./ccVars(j_o) )/sum(1./ccVars(j_o)) ) ...
        / segmentDuration;

      data_e(k,1) = ...
        sum( (inputData(j_e,k)./ccVars(j_e) )/sum(1./ccVars(j_e)) ) ...
        / segmentDuration;
    end
  end

  % sens ints ------------------------------------------------------------
  if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
    data_o = transpose(sum(inputData(j_o,:),1))*(segmentDuration^2);
    data_e = transpose(sum(inputData(j_e,:),1))*(segmentDuration^2);
  end                                                                            
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
  if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
    dataI  = inputData(1:end-1,:) * (segmentDuration^2);
    dataJ  = inputData(2:end,:) * (segmentDuration^2);
    dataIJ = 0.5*(woverlap4bar/w4bar)*0.5* ...
             transpose(sum(dataI,1)+sum(dataJ,1)); 
  end
                                                      
  % construct optimal weighting factors from the covariance matrix
  lambda_o = (C_ee - C_oe)/detC;
  lambda_e = (C_oo - C_oe)/detC;
  
  % cc stats --------------------------------------------------------------
  if ( strncmp(combineWhat, 'stats', length(combineWhat)) )
    combinedData = (combinedData_o*lambda_o + combinedData_e*lambda_e) ...
                   / (lambda_o + lambda_e) ;
  end

  % cc spectra ------------------------------------------------------------
  if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
    data = (data_o*lambda_o + data_e*lambda_e) / (lambda_o + lambda_e) ;
  end

  % sens ints -------------------------------------------------------------
  if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
    data = (data_o + data_e - 2*dataIJ)*((C_oo*C_ee)/detC);
  end

  % error bar -------------------------------------------------------------
  errorBar = sqrt(1/(lambda_o + lambda_e));

end

% construct freq series if necessary
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) | ...
     strncmp(combineWhat, 'sensints', length(combineWhat)) )
  combinedData = constructFreqSeries(data, flow, deltaF);
end

return

