function [combinedData, errorBar, numSegmentsTotal, inputData] = ...
  combineResultsFromFile(combineWhat, filename, segmentDuration, ...
                         doOverlap, window1, window2, ...
                         badGPSTimes, displayresults)
%
%  [combinedData, errorBar, numSegmentsTotal, inputData] = ...
%   combineResultsFromFile(combineWhat, filename, segmentDuration, ...
%                          doOverlap, window1, window2, ...
%                          badGPSTimes, displayresults)
%
%  read in cc stats, cc spectra or sensitivity integrands from a file,
%  then call combineResults.m to optimally combine them.
%  
%  Input:
%
%    combineWhat = 'stats', 'cmplxstats', 'spectra', 'sensints'
%
%    filename = name of file with cc stats, cc spectra, or sens ints
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    badGPSTimes = column vector with GPS times to veto
%    displayresults = 0 don't write results to screen (default = 1)
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
%    numSegmentsTotal = total number of data segments combined
%
%    inputData = data actually combined
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: combineResultsFromFile.m,v 1.17 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
if ( strncmp(combineWhat, 'stats', length(combineWhat)) ...
     || strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
  combinedData = 0;
else
  combinedData = constructFreqSeries(0,0,0);
end
errorBar = 0;
numSegmentsTotal = 0;

% initialize freqs (can be anything for cc stats)
flow = 0;
fhigh = 0;
deltaF = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in formatted data from a file, extracting relevant input data

% cc stats --------------------------------------------------------------
if ( strncmp(combineWhat, 'stats', length(combineWhat)) )
  [gpsTimes, stats, sigmas] = ...
    ctextread(filename, '%f%f%f\n', -1, 'commentstyle', 'matlab');
  if isempty(gpsTimes)
    inputData = 0;
    return;
  end

  ccSigmas = sigmas;
  inputData = stats;
end

% complex cc stats
if ( strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
  [gpsTimes, stats_real, stats_imag, sigmas] = ...
    ctextread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
  if isempty(gpsTimes)
    inputData = 0;
    return;
  end

  ccSigmas = sigmas;
  inputData = stats_real + 1i*stats_imag;
end

% cc spectra ------------------------------------------------------------
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
  [gpsTimes, sigmas, freqs, spectra_real, spectra_imag] = ...
    ctextread(filename, '%f%f%f%f%f\n', -1, 'commentstyle', 'matlab');
  if isempty(gpsTimes)
    inputData = [];
    return;
  end

  [gpsTimes, ind] = unique(gpsTimes);
  freqs        = unique(freqs);

  numFreqs = length(freqs);
  flow = freqs(1);
  fhigh = freqs(end);
  deltaF = freqs(2)-freqs(1);

  ccSigmas = sigmas(ind);
  inputData = reshape(spectra_real+1i*spectra_imag, numFreqs, length(gpsTimes));
  inputData = transpose(inputData);
end

% sens ints ------------------------------------------------------------
if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
  [gpsTimes, sigmas, freqs, integrands] = ...
    ctextread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
  if isempty(gpsTimes)
    inputData = [];
    return;
  end

  [gpsTimes, ind] = unique(gpsTimes);
  freqs        = unique(freqs);

  numFreqs = length(freqs);
  flow = freqs(1);
  fhigh = freqs(end);
  deltaF = freqs(2)-freqs(1);

  ccSigmas = sigmas(ind);
  inputData = reshape(integrands, numFreqs, length(gpsTimes));
  inputData = transpose(inputData);
end

% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
inputData = inputData(ind,:);
ccSigmas  = ccSigmas(ind);

% ignore bad gps times
[gpsTimes, ind] = setdiff(gpsTimes,badGPSTimes);
inputData = inputData(ind,:);
ccSigmas  = ccSigmas(ind);

% check that there's still exists data to combine
if isempty(gpsTimes) return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine the data
numSegmentsTotal = length(gpsTimes);
ii=1;
kk=1;
if ( strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
  combineWhatArg = 'stats';
else
  combineWhatArg = combineWhat;
end;

% Eliminate warnings in compiler
inputData_ovl = [];
ccSigmas_ovl = [];

if doOverlap==1

  while ii <= numSegmentsTotal

    clear inputData_ovl;
    clear ccSigmas_ovl;

    jj=1;
    inputData_ovl(jj,:) = inputData(ii,:);
    ccSigmas_ovl(jj,1)  = ccSigmas(ii);

    if ii~=numSegmentsTotal
      while gpsTimes(ii+1) == gpsTimes(ii)+segmentDuration/2
        ii=ii+1;
        jj=jj+1;
        inputData_ovl(jj,:) = inputData(ii,:);
        ccSigmas_ovl(jj,1)  = ccSigmas(ii);
        if ii==numSegmentsTotal
          break;
        end  
      end
    end

    [combinedData, errorBar, numSegments] = ...
      combineResults(combineWhatArg, inputData_ovl, ccSigmas_ovl, ...
                     segmentDuration, flow, fhigh, deltaF, ...
                     doOverlap, window1, window2);

    % apply appropriate factors of segmentDuration

    % error bar --------------------------------------------------------
    ccSigmas_nonovl(kk,1)  = errorBar*segmentDuration;
    % stats ------------------------------------------------------------
    if ( strncmp(combineWhat, 'stats', length(combineWhat)) ...
	 || strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
      inputData_nonovl(kk,:) = combinedData * segmentDuration;
    end
    % spectra ----------------------------------------------------------
    if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
      inputData_nonovl(kk,:) = transpose(combinedData.data) * segmentDuration;
    end
    % sens ints --------------------------------------------------------
    if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
      inputData_nonovl(kk,:) = transpose(combinedData.data) ...
                               / (segmentDuration^2);
    end

    ii=ii+1;
    kk=kk+1;

  end
   
  % combine cc stats, etc. from non-overlapping segments
  [combinedData, errorBar, numSegments] = ...
    combineResults(combineWhatArg, inputData_nonovl, ccSigmas_nonovl, ...
                   segmentDuration, flow, fhigh, deltaF, 0, 0, 0);

else
  [combinedData, errorBar, numSegmentsTotal] = ...
    combineResults(combineWhatArg, inputData, ccSigmas, segmentDuration, ...
                   flow, fhigh, deltaF, 0, 0, 0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write results to the screen
try 
  displayresults;
catch
  displayresults = 1;
end
if displayresults==0
  return
end

% cc stats --------------------------------------------------------------
if ( strncmp(combineWhat, 'stats', length(combineWhat)) )
  fprintf('Pt estimate = %e\n', combinedData);
end;

% complex cc stats ------------------------------------------------------
if ( strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
  if (imag(combinedData)>=0)
    fprintf('Pt estimate = %e + %e i\n',...
	    real(combinedData), imag(combinedData));
  else
    fprintf('Pt estimate = %e - %e i\n',...
	    real(combinedData), -imag(combinedData));
  end;
end;

% cc spectra ------------------------------------------------------------
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
  fprintf('Pt estimate = %e\n', ...
          2*sum(real(combinedData.data))*combinedData.deltaF);
end

% sens ints -------------------------------------------------------------
if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
  errorBarAlt = sqrt(1/sum(combinedData.data*combinedData.deltaF));
  fprintf('Alternate error bar calculation = %e\n', errorBarAlt);
end

fprintf('Error bar   = %e\n', errorBar);
fprintf('Number segs = %d\n', numSegmentsTotal);

return
