function [combinedData, errorBar, numSegmentsTotal, allData] = ...
  combineResultsFromJobs(combineWhat, filePrefix, fileSuffix, numJobs, ...
                         segmentDuration, doOverlap, window1, window2, ...
                         badGPSTimes, outputFileNamePrefix, ...
                         displayresults, figureNumber, figureLegend)
%
%  WARNING: This function will be deprecated!
%
%  [combinedData, errorBar, numSegmentsTotal, allData] = ...
%   combineResultsFromJobs(combineWhat, filePrefix, fileSuffix, numJobs, ...
%                          segmentDuration, doOverlap, window1, window2, ...
%                          badGPSTimes, outputFileNamePrefix, ...
%                          displayresults, figureNumber, figureLegend)
%
%  optimally combines cc stats, cc spectra, or sensitivity integrands
%  for a collection of jobs, by calling combineResultsFromFile repeatedly.
%
%  Input:
%
%    combineWhat = 'stats', 'spectra', 'sensints'
%
%    filePrefix, fileSuffix: prefix, suffix of filenames containing
%    the cc stat, cc spectra, or sensitivity integrand data
%    (note: the filename is assumed to be of the form
%    filePrefixNfileSuffix, where N runs from 1 to numJobs
%
%    numJobs = total number of jobs 
%
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    badGPSTimes = column vector with GPS times to veto
%    outputFileNamePrefix = prefix of filename that will contain
%      the integrand of the point estimate or sensitivity integrand;
%      prefix of filenames for figures if displaying results.
%    displayresults = 0 don't write results to screen (default = 1)
%    figureNumber = number of first figure (if displayresults = 1)
%    figureLegend = text for figure legend (if displayresults = 1)
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
%    allData = data actually combined
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: combineResultsFromJobs.m,v 1.17 2006-04-14 21:05:42 whelan Exp $
%
%  WARNING: This function will be deprecated!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddmmyyyyhhmmss  = datestr(now);

% default return values
if ( strncmp(combineWhat, 'stats', length(combineWhat)) ...
     || strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
else
  combinedData = constructFreqSeries(0,0,0);
end
errorBar = 0;
numSegmentsTotal = 0;

% initialize freqs (can be anything for cc stats)
flow = 0;
fhigh = 0;
deltaF = 0;

% initialize arrays
if ( strncmp(combineWhat, 'stats', length(combineWhat)) ...
     || strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
  combinedDataArray=zeros(numJobs, 1);
else
  combinedDataArray=struct('data',{},'flow',{},'deltaF',{},'symmetry',{});
end
errorBars   = zeros(numJobs,1);
numSegments = zeros(numJobs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over number of jobs

% WARNING "alldata" does not work on spectra yet;
allData = [];

for k = 1:numJobs
 
  fprintf('Analysing job %d\n', k);

  filename = [filePrefix num2str(k) fileSuffix];

  % check for missing job files
  if exist(filename)
    [combinedDataArray(k), errorBars(k), numSegments(k), dataFromFile]=...
      combineResultsFromFile(combineWhat, filename, segmentDuration, ...
                             doOverlap, window1, window2, badGPSTimes, 0);

%    allData = [allData; dataFromFile];

    if numSegments(k)==0
      fprintf('Job %d does not contain any cc spectra\n', k);
    end

  else
    fprintf('Missing job %d (filename %s)\n', k, filename)
  
  end

end

% eliminate zero entries associated with missing jobs, ...
ind=find(numSegments~=0);
combinedDataArray = combinedDataArray(ind);
errorBars = errorBars(ind);
numSegments = numSegments(ind);
numSegmentsTotal = sum(numSegments);
if numSegmentsTotal==0 return; end

% construct input data for final call to combineResults
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) | ...
     strncmp(combineWhat, 'sensints', length(combineWhat)) )

  numFreqs = length(combinedDataArray(1).data);
  deltaF = combinedDataArray(1).deltaF;
  flow  = combinedDataArray(1).flow;
  fhigh = flow + deltaF*(numFreqs-1);
  freqs = flow + deltaF*transpose([0:numFreqs-1]);  
                                                                              
  % extract data from the freq series structure
  inputData = zeros(length(numSegments), numFreqs);
  for k=1:length(numSegments)

    % different normalization for spectra and sens ints
    if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
      inputData(k,:)=transpose(combinedDataArray(k).data)*segmentDuration;
    else
      inputData(k,:)=transpose(combinedDataArray(k).data)/(segmentDuration^2);
    end

  end

else
  
  % cc stats
  inputData = combinedDataArray*segmentDuration;

end    

% error bars
ccSigmas = errorBars * segmentDuration;

% calculate final results (no overlap between segments from diff jobs)
if ( strncmp(combineWhat, 'cmplxstats', length(combineWhat)) )
[combinedData, errorBar, ignore] = ...
  combineResults('stats', inputData, ccSigmas, segmentDuration, ... 
                 flow, fhigh, deltaF, 0, 0, 0);
else
[combinedData, errorBar, ignore] = ...
  combineResults(combineWhat, inputData, ccSigmas, segmentDuration, ... 
                 flow, fhigh, deltaF, 0, 0, 0);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write combine cc spectra, sens ints to file 

% cc spectra ------------------------------------------------------------
if ( strncmp(combineWhat, 'spectra', length(combineWhat)) )
  filename = [outputFileNamePrefix '_ptEstIntegrand.dat'];
  fid = fopen(filename, 'w');
  fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
  fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
  fprintf(fid, '%s\t%s\t%s\t%s\n', ...
          '%error bar', 'freq (Hz)', ...
          'pt est int (re)', 'pt est int (im)');
  for ii=1:numFreqs
    fprintf(fid, '%e\t%e\t%e\t%e\n', errorBar, freqs(ii), ...
            real(combinedData.data(ii)), imag(combinedData.data(ii)));
  end 
  fclose(fid);
end

% sens ints -------------------------------------------------------------
if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )
  filename = [outputFileNamePrefix '_sensIntegrand.dat'];
  fid = fopen(filename, 'w');
  fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
  fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
  fprintf(fid, '%s\t%s\t%s\n', ...
          '%error bar', 'freq (Hz)', ...
          'sens integrand');
  for ii=1:numFreqs
    fprintf(fid, '%e\t%e\t%e\n', errorBar, freqs(ii), combinedData.data(ii));
  end 
  fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display final results if desired
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

end

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

  % plots
  figure(figureNumber)
  plot(freqs, real(combinedData.data));
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Real(Integrand of point estimate)', 'FontSize', 14);
  legend(figureLegend);
  filename = [outputFileNamePrefix '_ptEstIntegrand_real'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');

  figure(figureNumber+1)
  plot(freqs, imag(combinedData.data));
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Imag(Integrand of point estimate)', 'FontSize', 14);
  legend(figureLegend);
  filename = [outputFileNamePrefix '_ptEstIntegrand_imag'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');
                                                                              
  figure(figureNumber+2)
  plot(freqs, abs(combinedData.data));
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Abs(Integrand of point estimate)', 'FontSize', 14);
  legend(figureLegend);
  filename = [outputFileNamePrefix '_ptEstIntegrand_abs'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');
  
  figure(figureNumber+3)
  cumpointest = cumsum(2*real(combinedData.data)*deltaF);
  plot(freqs, cumpointest);
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Cumulative point estimate', 'FontSize', 14);
  legend(figureLegend);
  filename = [outputFileNamePrefix '_ptEstIntegrand_cum'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');

end

% sens ints -------------------------------------------------------------
if ( strncmp(combineWhat, 'sensints', length(combineWhat)) )

  errorBarAlt = sqrt(1/sum(combinedData.data*combinedData.deltaF));
  fprintf('Alternate error bar calculation = %e\n', errorBarAlt);

  % plots
  figure(figureNumber)
  plot(freqs, combinedData.data);
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Sensitivity Integrand', 'FontSize', 14);
  legend(figureLegend);
  filename = [outputFileNamePrefix '_sensIntegrand'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');
  
  figure(figureNumber+1)
  cumsens = cumsum(combinedData.data);
  cumsens = cumsens/cumsens(end);
  plot(freqs, cumsens);
  grid on;
  xlabel('f (Hz)', 'FontSize', 14);
  ylabel('Cumulative sensitivity', 'FontSize', 14);
  legend(figureLegend, 2);
  filename = [outputFileNamePrefix '_sensIntegrand_cum'];
  print(gcf, '-depsc2', filename);
  saveas(gcf, filename, 'fig');
  
end

fprintf('Error bar   = %e\n', errorBar);
fprintf('Number segs = %d\n', numSegmentsTotal);

return
