function [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
          numSegmentsTotal] = ...
  combineResultsFromMultipleJobs(ccStatsFilePrefix, ccSpectraFilePrefix, ...
                                 sensIntsFilePrefix, fileSuffix, ...
                                 numJobs, segmentDuration, badGPSTimes, ...
                                 doRenormalize, modifyFilter, ...
                                 doOverlap, window1, window2, ...
                                 outputFileNamePrefix, displayResults, ...
                                 figureNumber, figureLegend, complexflag, ...
                                 recoveryFile)
%
%  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
%   numSegmentsTotal] = ...
%  combineResultsFromMultipleJobs(ccStatsFilePrefix, ccSpectraFilePrefix, ...
%                                 sensIntsFilePrefix, fileSuffix, ...
%                                 numJobs, segmentDuration, badGPSTimes ...
%                                 doRenormalize, modifyFilter, ...
%                                 doOverlap, window1, window2, ...
%                                 outputFileNamePrefix, displayResults, ... 
%                                 figureNumber, figureLegend, complexflag, ...
%                                 recoveryFile)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands for a collection of jobs, by calling 
%  combineResultsFromSingleJob repeatedly.
%
%  Input:
%
%    ccStatsFilePrefix, ccSpectraFilePrefix, sensIntsFilePrefix: 
%    prefixes of filenames containing the cc stats, cc spectra, and 
%    sensitivity integrand data
%
%    fileSuffix: suffix of filenames for the above data
%    NOTE: the filename is assumed to be of the form 
%    filePrefixNfileSuffix, where N runs from 1 to numJobs
%
%    numJobs = total number of jobs
%
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    badGPSTimes = column vector with GPS times to veto
%
%    doRenormalize = 0 don't renormalize for a different optimal filter
%                  = 1 renormalize for a different optimal filter
%    modifyFilter = row vector containing values for modifying the optimal 
%      filter
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%
%    outputFileNamePrefix = prefix of filename that will contain
%      the integrand of the point estimate and sensitivity integrand;
%      prefix of filenames for figures if displaying results.
%    displayResults = 0 don't write results to screen (default = 1)
%    figureNumber = number of first figure (if displayResults = 1)
%    figureLegend = text for figure legend (if displayResults = 1)
%
%    complexflag = 0 means CC stats real; 1 means complex
%
%    if recoveryFile is given, combineResultsFromMultipleJobs will resume
%    from the last iteration saved.
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
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk or john.whelan@ligo.org
%
%  $Id: combineResultsFromMultipleJobs.m,v 1.12 2006-09-11 23:16:42 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyyymmddhhMMss  = datestr(now,31);
yyyymmdd = datestr(now,29);

% check complex flag
try
  complexflag;
catch
  complexflag = false;
end;


% Check whether or not to save/load recovery file
try
  recoveryFile;
catch
  recoveryFile = '';
end;

% default return values
ptEstimate = 0;
errorBar = 0;
combinedPtEstInt = constructFreqSeries(0,0,0);
combinedSensInt = constructFreqSeries(0,0,0);
numSegmentsTotal = 0;

% initialize arrays
ptEstimates = zeros(numJobs,1);
errorBars = zeros(numJobs,1);
combinedPtEstInts = struct('data',{},'flow',{},'deltaF',{},'symmetry',{});
combinedSensInts = struct('data',{},'flow',{},'deltaF',{},'symmetry',{});
numSegments = zeros(numJobs,1);
ccStatsRenormalized = [];
ccSigmasRenormalized = [];
goodGPSTimes = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strncmp(ccSpectraFilePrefix,'<none>', length(ccSpectraFilePrefix)))
  calculateSpectra = false;
else
  calculateSpectra = true;
end;

startpos = 1;
% If there's already saved work, load it and continue; numJobs may have been
% updated, though, so take the new value.
if (~strcmp(recoveryFile,''))
  newnumJobs = numJobs;
  try
    load(recoveryFile);
    fprintf(1,['Resuming state from ' recoveryFile '.\n']);
    startpos = k;
  catch
    fprintf(1,['Unable to load recovery file' recoveryFile...
               '.  Will checkpoint.\n']);
  end;
  numJobs = newnumJobs;
  
  clear newnumJobs % Careful not to re-save this value later
end;

for k = startpos:numJobs
  % If recoveryFile is given, save state every 50th iteration
  if (mod(k,50)==0 && ~strcmp(recoveryFile,''))
    save(recoveryFile)
    fprintf(1,['Checkpointed through job ' num2str(k) '.\n'])
  end;

  statsfilename = [ccStatsFilePrefix num2str(k) fileSuffix];
  spectrafilename = [ccSpectraFilePrefix num2str(k) fileSuffix];
  sensintsfilename = [sensIntsFilePrefix num2str(k) fileSuffix];

  if (~exist(statsfilename))
    fprintf('Missing job %d (filename %s)\n', k, statsfilename)
    continue;
    %return;
  end;

  if (calculateSpectra && ~exist(spectrafilename))
    fprintf('Missing job %d (filename %s)\n', k, spectrafilename)
    continue;
    %return;
  end;

  if (~exist(sensintsfilename))
    fprintf('Missing job %d (filename %s)\n', k, sensintsfilename)
    continue;
    %return;
  end;

  % read in data from files
  [ccStats, ccSigmas, totalGPSTimes] = ...
    readCCstatsCCsigmasFromFile(statsfilename, complexflag);
  [sensInts,  flow, deltaF, totalGPSTimes] = ...
    readSensIntsFromFile(sensintsfilename);

  % nvf: If the last job is empty, sometimes vars get left in a bad state.
  % Save it here in good loop.  Pull it up at the end of iteration if bad loop.
  if length(flow)~=0
    flowold = flow;
    deltaFold = deltaF;
  end;

  if calculateSpectra
    [ccSpectra, flow, deltaF, totalGPSTimes] = ...
        readCCspectraFromFile(spectrafilename);
  else
    ccSpectra=sensInts; % dummy value to allow future calculations
  end;
  if length(totalGPSTimes)~=0

    % construct fhigh
    fhigh = flow + deltaF*(length(ccSpectra(1,:))-1);

    % ignore bad gps times
    [ccStats, gpsTimes] = ...
      ignoreBadGPSTimes(ccStats, totalGPSTimes, badGPSTimes);
    [ccSigmas, gpsTimes] = ...
      ignoreBadGPSTimes(ccSigmas, totalGPSTimes, badGPSTimes);
    [ccSpectra, gpsTimes] = ...
        ignoreBadGPSTimes(ccSpectra, totalGPSTimes, badGPSTimes);
    [sensInts,  gpsTimes] = ...
    ignoreBadGPSTimes(sensInts, totalGPSTimes, badGPSTimes);

    % renormalize data if desired
    if doRenormalize
      fprintf('renormalizing the data\n');
      for jj = 1:length(gpsTimes)

        [ccStats(jj,1), ccSigmas(jj,1), ccSpectra(jj,:), sensInts(jj,:)] = ...
          renormalizeData(ccSpectra(jj,:), sensInts(jj,:), modifyFilter, ...
                            flow, fhigh, deltaF);

      end
    end

    % construct arrays of renormalised cc stats and theor sigmas for 
    % the good GPS times
    ccStatsRenormalized = [ccStatsRenormalized; ccStats(:,1)];
    ccSigmasRenormalized = [ccSigmasRenormalized; ccSigmas(:,1)];
    goodGPSTimes = [goodGPSTimes; gpsTimes];

    % combine renormalized data
    [ptEstimates(k), errorBars(k), combinedPtEstInts(k), ...
      combinedSensInts(k), numSegments(k)] = ...
      combineResultsFromSingleJob(ccStats, ccSigmas, ccSpectra, sensInts, ...
                                  gpsTimes, segmentDuration, ...
                                  flow, fhigh, deltaF, ...
                                  doOverlap, window1, window2, 0);
  else

    numSegments(k)=0;

  end

  if numSegments(k)==0
    fprintf('Job %d does not contain any data\n', k);
    flow = flowold;
    deltaF = deltaFold;
  end

end % loop over numJobs

if (~strcmp(recoveryFile,''))
  k=k+1;
  save(recoveryFile)
  fprintf(1,['Checkpointed through job ' num2str(k) '.\n'])
end;

% eliminate zero entries associated with missing jobs, ...
ind=find(numSegments~=0);
ptEstimates = ptEstimates(ind);
errorBars = errorBars(ind);
combinedPtEstInts = combinedPtEstInts(ind);
combinedSensInts = combinedSensInts(ind);
numSegments = numSegments(ind);

numSegmentsTotal = sum(numSegments);
if numSegmentsTotal==0 return; end

% extract data from the freq series structures
numFreqs = length(combinedPtEstInts(1).data);
freqs = flow + deltaF*transpose([0:numFreqs-1]);
ccSpectra = zeros(length(numSegments), numFreqs);
sensInts = zeros(length(numSegments), numFreqs);
for k=1:length(numSegments)
  ccSpectra(k,:)=transpose(combinedPtEstInts(k).data)*segmentDuration;
  sensInts(k,:)=transpose(combinedSensInts(k).data)/(segmentDuration^2);
end    

ccStats = ptEstimates*segmentDuration;
ccSigmas = errorBars*segmentDuration;

% calculate final results (no overlap between segments from diff jobs)
[ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ignore] = ...
  combineResultsNew(ccStats, ccSigmas, ccSpectra, sensInts, ...
                    segmentDuration, flow, fhigh, deltaF, 0, 0, 0, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write renormalised cc stats, theor sigmas for good GPS times to file

filename = [outputFileNamePrefix '_ccstats_renormalized.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', yyyymmddhhMMss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
if complexflag
  fprintf(fid, '%s\t%s\t%s\t%s\n', ...
          '%start sec', 'Re(CC stat)', 'Im(CC stat)', 'theor sigma');
  for ii=1:length(goodGPSTimes)
    fprintf(fid, '%d\t%e\t%e\t%e\n', ...
	    goodGPSTimes(ii), real(ccStatsRenormalized(ii)), ...
	    imag(ccStatsRenormalized(ii)), ccSigmasRenormalized(ii));
  end;
else
  fprintf(fid, '%s\t%s\t%s\n', ...
          '%start sec', 'CC statistic', 'theor sigma');
  for ii=1:length(goodGPSTimes)
    fprintf(fid, '%d\t%e\t%e\n', ...
	    goodGPSTimes(ii), ccStatsRenormalized(ii), ...
	    ccSigmasRenormalized(ii));
  end;
end
fclose(fid);
                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write combined cc spectrum, sens int to file

% cc spectrum
filename = [outputFileNamePrefix '_ptEstIntegrand.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', yyyymmddhhMMss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\t%s\t%s\n', ...
        '%error bar', 'freq (Hz)', 'pt est int (re)', 'pt est int (im)');
for ii=1:numFreqs
  fprintf(fid, '%e\t%e\t%e\t%e\n', errorBar, freqs(ii), ...
          real(combinedPtEstInt.data(ii)), imag(combinedPtEstInt.data(ii)));
end
fclose(fid);
                                                                                
% sens int
filename = [outputFileNamePrefix '_sensIntegrand.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', yyyymmddhhMMss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\t%s\n', ...
        '%error bar', 'freq (Hz)', 'sens integrand');
for ii=1:numFreqs
  fprintf(fid, '%e\t%e\t%e\n', errorBar, freqs(ii), combinedSensInt.data(ii));
end
fclose(fid);

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

% Check whether we're running compiled matlab
try
  deployedFlag = isdeployed;
catch
  deployedFlag = 0;
end;

if complexflag
  ptEstimateAlt = sum((combinedPtEstInt.data)*combinedPtEstInt.deltaF);
else
  ptEstimateAlt = 2*sum(real(combinedPtEstInt.data)*combinedPtEstInt.deltaF);
end;
errorBarAlt = sqrt(1/sum(combinedSensInt.data*combinedSensInt.deltaF));

%fprintf('combineResultsFromMultipleJobs.m results -----------------------\n');
if complexflag
  fprintf('Pt estimate = %e + (%e) i\n', real(ptEstimate), imag(ptEstimate));
else
  fprintf('Pt estimate = %e\n', ptEstimate);
end;
fprintf('Error bar   = %e\n', errorBar);
if complexflag
  fprintf('Alternate point estimate calculation = %e + (%e) i\n',...
	  real(ptEstimateAlt), imag(ptEstimateAlt));
else
  fprintf('Alternate point estimate calculation = %e\n', ptEstimateAlt);
end;
fprintf('Alternate error bar calculation = %e\n', errorBarAlt);
fprintf('Number segs = %d\n', numSegmentsTotal);
%fprintf('combineResultsFromMultipleJobs.m results -----------------------\n');
                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cc spectrum
figure(figureNumber)
plot(freqs, real(combinedPtEstInt.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Real(Integrand of point estimate)', 'FontSize', 14);
legend(figureLegend);
title(['ptEstIntegrand\_real ' yyyymmdd]);
filename = [outputFileNamePrefix '_ptEstIntegrand_real'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;
                                                                                
figure(figureNumber+1)
plot(freqs, imag(combinedPtEstInt.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Imag(Integrand of point estimate)', 'FontSize', 14);
legend(figureLegend);
title(['ptEstIntegrand\_imag ' yyyymmdd]);
filename = [outputFileNamePrefix '_ptEstIntegrand_imag'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;

figure(figureNumber+2)
plot(freqs, abs(combinedPtEstInt.data));
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Abs(Integrand of point estimate)', 'FontSize', 14);
legend(figureLegend);
title(['ptEstIntegrand\_abs ' yyyymmdd]);
filename = [outputFileNamePrefix '_ptEstIntegrand_abs'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;
                                                                                
figure(figureNumber+3)
cumpointest = cumsum(2*real(combinedPtEstInt.data)*deltaF);
plot(freqs, cumpointest);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative point estimate', 'FontSize', 14);
legend(figureLegend);
title(['ptEstIntegrand\_cum ' yyyymmdd]);
filename = [outputFileNamePrefix '_ptEstIntegrand_cum'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;
                                                                                
% sens int
figure(figureNumber+4)
plot(freqs, combinedSensInt.data);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Sensitivity Integrand', 'FontSize', 14);
legend(figureLegend);
title(['sensIntegrand ' yyyymmdd]);
filename = [outputFileNamePrefix '_sensIntegrand'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;
                                                                                
figure(figureNumber+5)
cumsens = cumsum(combinedSensInt.data);
cumsens = cumsens/cumsens(end);
plot(freqs, cumsens);
grid on;
xlabel('f (Hz)', 'FontSize', 14);
ylabel('Cumulative sensitivity', 'FontSize', 14);
legend(figureLegend, 2);
title(['sensIntegrand\_cum ' yyyymmdd]);
filename = [outputFileNamePrefix '_sensIntegrand_cum'];
print(gcf, '-depsc2', filename);
print(gcf, '-dpng', filename);
saveas(gcf, filename, 'fig');
if deployedFlag
    delete(gcf);
end;
                                                                                
return
