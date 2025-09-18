function [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
          numSegmentsTotal] = ...
  combineResultsFromMultipleJobs_H1H2(commonFilePrefix, fileSuffix, ...
                                 StarT, EnD, segmentDuration, badGPSTimes, ...
                                 doRenormalize, modifyFilter, ...
                                 doOverlap, window1, window2, ...
                                 outputFileNamePrefix, displayResults, ...
                                 figureNumber, figureLegend, ...
                                 runFlag, jobsFile, complexflag)
%
%  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
%   numSegmentsTotal] = ...
%  combineResultsFromMultipleJobs(commonFilePrefix, fileSuffix, ...
%                                 startJob, endJob, segmentDuration, badGPSTimes ...
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
%    startJob = Number of starting Job
%
%    endJob = Number of end Job
%
%    segmentDuration = length of analysis segment in sec (typically 60)
%
%    badGPSTimes = column vector with GPS times to veto
%
%    doRenormalize = 0 don't renormalize for a different optimal filter or freq mask
%                  = 1 renormalize for a different optimal filter or freq mask
%    modifyFilter = row vector containing values for modifying the optimal 
%                   filter or row vector of zeros and ones for freq mask
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
% $Id:combineResultsFromMultipleJobs.m,v 1.15 2010/03/05 23:16:42 shivaraj Exp$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyyymmddhhMMss  = datestr(now,31);
yyyymmdd = datestr(now,29);

% check complex flag
if(exist('complexflag','var')==0)
  complexflag = false;
end;

% Flag to indiate whether the given Start and EnD are GPS times are Job
% numbers; 0 - Job numbers, 1 - GPS times
if(exist('runFlag','var')==0)
  runFlag = 0;
end;

% File prefix for ccstat, spectra and sensint
ccStatsFilePrefix = [commonFilePrefix '_ccstats.job'];
ccSpectraFilePrefix = [commonFilePrefix '_ccspectra.job'];
sensIntsFilePrefix = [commonFilePrefix '_sensints.job'];

% default return values
ptEstimate = 0;
errorBar = 0;
combinedPtEstInt = constructFreqSeries(0,0,0);
combinedSensInt = constructFreqSeries(0,0,0);
numSegmentsTotal = 0;

% Determining the number of Jobs that will be enough to contain the given time
% interval. By this method we can avoid running over all the jobs.
if runFlag
  Job_details=load(jobsFile);
  AllStartTimes = Job_details(:,2);
  AllEndTimes = Job_details(:,3);
  cut = AllEndTimes-segmentDuration>=StarT & ...
         AllStartTimes+segmentDuration<=EnD;
  jobno_series = [1:length(AllStartTimes)];
  Probable_jobnos = jobno_series(cut);%this estimate is approx, but conservative
  if(length(Probable_jobnos)~=0)
    startJob = Probable_jobnos(1);
    endJob = Probable_jobnos(end);
  else
   fprintf('There is no job in the given interval\n');
   return;
  end
elseif ~runFlag
 startJob = StarT;
 endJob = EnD;
end

% initialize arrays
numJobs = endJob - startJob + 1;
ptEstimates = zeros(numJobs,1);
errorBars = zeros(numJobs,1);
combinedPtEstInts = struct('data',{},'flow',{},'deltaF',{},'symmetry',{});
combinedSensInts = struct('data',{},'flow',{},'deltaF',{},'symmetry',{});
numSegments = zeros(numJobs,1);
ccStatsRenormalized = []; %for additional freq mask or filter
ccSigmasRenormalized = []; %for additional freq mask or filter
goodGPSTimes = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = startJob:endJob
  fprintf('Analysing Job no %d \n',k);
  statsfilename = [ccStatsFilePrefix num2str(k) fileSuffix];
  spectrafilename = [ccSpectraFilePrefix num2str(k) fileSuffix];
  sensintsfilename = [sensIntsFilePrefix num2str(k) fileSuffix];

  if (~exist(statsfilename))
    fprintf('Missing job %d (filename %s)\n', k, statsfilename)
    continue;
    %return;
  end;

  if (~exist(spectrafilename))
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
  [sensInts, flow, deltaF, totalGPSTimes] = ...
          readSensIntsFromFile(sensintsfilename);
  [ccSpectra, flow, deltaF, totalGPSTimes] = ...
        readCCspectraFromFile(spectrafilename);
  
  % If the last job is empty, then flow and deltaF becomes zero and it creates
  % problem. So we save a good value here and pull it up at the end when 
  % all results are combined.
  if length(flow)~=0
    flowgood = flow;
    deltaFgood = deltaF;
  end;

  if length(totalGPSTimes)~=0

    % construct fhigh
    fhigh = flow + deltaF*(length(ccSpectra(1,:))-1);

    % ignore bad gps times
    [gpsTimes, ind] = setdiff(totalGPSTimes, badGPSTimes);
    ccStats = ccStats(ind,:);
    ccSigmas = ccSigmas(ind,:);
    ccSpectra = ccSpectra(ind,:);
    sensInts = sensInts(ind,:);

    % Selecting only the segments that are in the given time interval; 
    % if GPS times are given 
    if(runFlag == 1)
      cut_gps = gpsTimes >= StarT & gpsTimes+segmentDuration <= EnD; 
      gpsTimes = gpsTimes(cut_gps,:);
      ccStats = ccStats(cut_gps,:);
      ccSigmas = ccSigmas(cut_gps,:);
      ccSpectra = ccSpectra(cut_gps,:);
      sensInts = sensInts(cut_gps,:);
    end

    % renormalize data if desired; with additional Freq Mask or filter
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

    % combine data
    [ptEstimates(k), errorBars(k), combinedPtEstInts(k), ...
      combinedSensInts(k), numSegments(k)] = ...
      combineResultsFromSingleJob(ccStats, ccSigmas, ccSpectra, sensInts, ...
                                  gpsTimes, segmentDuration, ...
                                  flow, fhigh, deltaF, ...
                                  doOverlap, window1, window2, 0);
  else

    numSegments(k)=0;
    fprintf('Job %d does not contain any data\n', k);

  end

end % loop over numJobs

% eliminate zero entries associated with missing jobs, ...
ind=numSegments>0;
ptEstimates = ptEstimates(ind);
errorBars = errorBars(ind);
combinedPtEstInts = combinedPtEstInts(ind);
combinedSensInts = combinedSensInts(ind);
numSegments = numSegments(ind);

numSegmentsTotal = sum(numSegments);
if numSegmentsTotal==0 
 fprintf('No data avilable or data is not enough \n'); 
 return; 
end

%In case the last job is empty this is needed
if length(flow)==0
    flow = flowgood;
    deltaF = deltaFgood;
end;

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

if doRenormalize
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
end                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write combined cc spectrum, sens int to file

if doRenormalize
 outputFileNamePrefix = [outputFileNamePrefix '_renormalized'];
end

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
if(exist('displayResults','var')==0)
  displayResults = 1;
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
if displayResults==0
  return
end

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
