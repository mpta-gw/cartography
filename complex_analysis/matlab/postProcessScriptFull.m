% postProcessScriptFull - script for post processing the data using the same
% inputs as stochastic.m, plus some additional parameters used only for
% post-processing.
%
% This script is intended to be run after first generating individual results
% for each job. Concatenation of individual job files is done automatically
% as part of this script using Unix calls to 'cat'. The complete set of steps
% needed to run the stochastic pipeline is:
%
% 1. Generate list of segments using SegWizard and save to the jobs file
% 2. Edit parameters file if necessary
% 3. Run stochastic_pipe.tclsh to create frame cache for all jobs
% 4. Run stochastic.m
% 5. Run postProcessScriptFull
%
% Inputs:
%   paramsFile - parameters file
%   jobsFile - list of science segments
%   outputDir - directory where post-processing results will go
%   dSigmaCut - max allowable relative difference between naive and theoretical
%               sigmas
%   largeSigmaCutoff
%   doRenormalize
%   modifyFilter
%   displayResults - set to "false" to suppress display of plots
%   badGPSTimesFile - file to which you wish to save badGPSTimes (optional)
%
% NOTE: when running the compiled postProcessScriptFull it is not possible to
% use a modified filter - modifyFilter must always be 0.
function postProcessScriptFull(paramsFile, jobsFile, outputDir, dSigmaCut, largeSigmaCutoff, doRenormalize, modifyFilter, displayResults, badGPSTimesFile)

  datestamp  = datestr(now,31);

  if (ischar(displayResults))
    displayResults = str2num(displayResults);
  end;

  if (ischar(dSigmaCut))
    dSigmaCut = str2num(dSigmaCut);
  end;

  if (ischar(largeSigmaCutoff))
    largeSigmaCutoff = str2num(largeSigmaCutoff);
  end;

  if (ischar(doRenormalize))
    doRenormalize = str2num(doRenormalize);
  end;

  if (ischar(modifyFilter))
    modifyFilter = str2num(modifyFilter);
  end;

  if (~displayResults)
    set(0, 'DefaultFigureVisible', 'off');
    set(0, 'DefaultAxesVisible', 'off');
  end;

  if (~exist(outputDir, 'dir'))
    s=sprintf('Output dir %s does not exist', outputDir);
    error(s);
  end;

  try
    badGPSTimesFile;
  catch
    badGPSTimesFile='';
  end

  sigmaMinRatio = 1 - dSigmaCut;
  sigmaMaxRatio = 1 + dSigmaCut;

  % These appear to be the same for all jobs
  fileSuffix = '.trial1.dat';
  figureNumber = 1;

  % read in params structure from a file
  params = readParamsFromFile(paramsFile);
  try
    doCombine = params.doCombine;
  catch
    doCombine = false;
  end;
  try  % Only set doBadGPSTimes if cuts have already been made in stochastic.m
    params.minDSigRatio;
    params.maxDSigRatio;
    doBadGPSTimes = params.doBadGPSTimes;
  catch
    doBadGPSTimes = false;
  end;
  try
    writeNaiveSigmasToFiles = params.writeNaiveSigmasToFiles;
  catch
    writeNaiveSigmasToFiles = false;
  end;
  try
    writeSpectraToFiles = params.writeSpectraToFiles;
  catch
    writeSpectraToFiles = false;
  end;
  try
    writeSensIntsToFiles = params.writeSensIntsToFiles;
  catch
    writeSensIntsToFiles = false;
  end;
  try
    doOverlap = params.doOverlap;
  catch
    doOverlap = false;
  end;
  try
    segmentDuration = params.segmentDuration;
  catch
    error('segmentDuration parameter not set');
  end;
  try
    resampleRate1 = params.resampleRate1;
  catch
    error('resampleRate1 parameter not set');
  end;
  try
    resampleRate2 = params.resampleRate2;
  catch
    error('resampleRate2 parameter not set');
  end;
  try
    outputFilePrefix = params.outputFilePrefix;
  catch
    error('outputFilePrefix parameter not set');
  end;
  try
    ifo1 = params.ifo1;
  catch
    error('ifo1 parameter not set');
  end;
  try
    ifo2 = params.ifo2;
  catch
    error('ifo2 parameter not set');
  end;
  clear params;

  % Prefix for output from post-processing
  outputFileNamePrefix = [ outputDir '/' ifo1 ifo2 ];
  figureLegend = [ ifo1 '-' ifo2 ];

  % read in job start time and duration from a file so we can count jobs and segments
  [ignore1, startTimes, ignore2, jobDurations] = ...
    textread(jobsFile, '%n %n %n %n', -1, 'commentstyle', 'matlab');
  numJobs = length(startTimes);
  numSegments = sum(jobDurations) / segmentDuration;
  clear ignore1 ignore2;

  naivesigmasfileprefix = [ outputFilePrefix '_naivesigmas.job' ];
  
  if doCombine
    statsfileprefix = [ outputFilePrefix '_combined_ccstats.job' ];
    spectrafileprefix = [ outputFilePrefix '_combined_ccspectra.job' ];
    sensintsfileprefix = [ outputFilePrefix '_combined_sensints.job' ];
  else
    statsfileprefix = [ outputFilePrefix '_ccstats.job' ];
    spectrafileprefix = [ outputFilePrefix '_ccspectra.job' ];
    sensintsfileprefix = [ outputFilePrefix '_sensints.job' ];
  end

  numPoints1 = segmentDuration*resampleRate1;
  numPoints2 = segmentDuration*resampleRate2;
  window1 = tukeywin(numPoints1, segmentDuration*resampleRate1);
  window2 = tukeywin(numPoints2, segmentDuration*resampleRate2);

  concatNaiveSigmas = [ outputFilePrefix '_concatNaiveSigmas.dat' ];
  concatCCStats = [ outputFilePrefix '_concatCCStats.dat' ];

% Print datestamp
unix([ 'echo "' datestamp '" > ' outputFilePrefix '_datestamp.txt']);

% Concatenate files from all jobs, one job at a time to prevent
% shell overflow. Note "find" can't be used easily because we don't
% have the root directory where output files are kept, just a prefix.
% We could parse the prefix but that would be tedious.
% Job 1 is done first to create the concat file, then the others are
% appended to it
job = 1;
tail = [ num2str(job) '.trial*.dat' ];
unix([ 'cat ' naivesigmasfileprefix tail ' > ' concatNaiveSigmas ]);
unix([ 'cat ' statsfileprefix tail ' > ' concatCCStats ]);
for job = 2:numJobs
  tail = [ num2str(job) '.trial*.dat' ];
  unix([ 'cat ' naivesigmasfileprefix tail ' >> ' concatNaiveSigmas ' 2> /dev/null']);
  unix([ 'cat ' statsfileprefix tail ' >> ' concatCCStats ' 2> /dev/null']);
end;

  if (writeNaiveSigmasToFiles)
    % relative sigma cuts
    [gpsTimes, sigmas, naiveSigmas, badGPSTimes_abssigmacut, ...
     goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
       compareSigmasFromFile(concatNaiveSigmas, sigmaMinRatio, sigmaMaxRatio);
      
    % absolute sigma cuts
    badGPSTimes_largesigmacut = largeSigmas(concatNaiveSigmas,largeSigmaCutoff);

    % union cuts
    badGPSTimes = [badGPSTimes_abssigmacut; badGPSTimes_largesigmacut];
    badGPSTimes = sort(unique(badGPSTimes));

    % report results
    if displayResults
      fprintf(1,[datestamp '\n']);
      fprintf(['Number of segments lost to abs sigma cut = %d,' ...
        ' fraction = %g\n'], ...
              length(badGPSTimes_abssigmacut), ...
              length(badGPSTimes_abssigmacut)/numSegments);
      fprintf(['Number of segments lost to large sigma cut = %d,' ...
               ' fraction = %g\n'], ...
              length(badGPSTimes_largesigmacut), ...
              length(badGPSTimes_largesigmacut)/numSegments);
      fprintf('Total number of segments lost to cuts = %d, fraction = %g\n', ...
              length(badGPSTimes), length(badGPSTimes)/numSegments);
      fprintf('Total time remaining after cuts = %5.1f days\n', ...
              (numSegments-length(badGPSTimes))*segmentDuration/(3600*24));
    end;
    clear badGPSTimes_abssigmacut badGPSTimes_largesigmacut;

    % save results 
    if length(badGPSTimesFile) > 0
      % in mat-file format if possible
      if strcmp(badGPSTimesFile(length(badGPSTimesFile)-3:end),'.mat')
        save(badGPSTimesFile,'badGPSTimes','-v6')
      else
        fid = fopen(badGPSTimesFile,'w');
        for ii = 1:length(badGPSTimes)
          fprintf(fid, '%d\n', badGPSTimes(ii));
        end;
        fclose(fid);
      end;
    else
      warning('No badGPSTimesFile file provided; not saving badGPSTimes.');
    end;
  else % if writeNaiveSigmasToFiles
    warning('Naive sigmas not written - sigma cuts will not be performed');
    badGPSTimes = [];
  end; % if writeNaiveSigmasToFiles

  % Ok, done with pre-post-processing.  Let's start the post-processing.
  % For keeping plots straight, figure #s are as numbered in Vuk's elog post at
  % http://ldas-sw.ligo.caltech.edu/ilog/pub/ilog.cgi?group=stochastic&date_to_view=08/31/2005&anchor_to_scroll_to=2005:08:31:12:10:03-vmandic

  % Figure 1 (+ two new figures)
  runningPointEstimate(concatCCStats, badGPSTimes, resampleRate1, ...
                       segmentDuration, outputFileNamePrefix, figureLegend, ...
                       doOverlap, displayResults);
  fprintf('** Done runningPointEstimate\n');

  % Figures 2-7
  [ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, numSegmentsTotal] = ...
  combineResultsFromMultipleJobs(statsfileprefix, spectrafileprefix, ...
          sensintsfileprefix, fileSuffix, numJobs, segmentDuration, ...
          badGPSTimes, doRenormalize, modifyFilter, doOverlap, window1, ...
          window2, outputFileNamePrefix, displayResults, figureNumber, ...
          figureLegend, 0, [outputFilePrefix '_resumeCRfMJ.mat']);
  fprintf('** Done combineResultsFromMultipleJobs\n');

  % Figures 9-17 (+ one new figure)
  % (there are only ten PanelPlots... "2" and "10" are missing)
  % doOverlap=1, DOFscalefactor = 1/(1+3/35) for 50% overlapping Hann
  StatisticalAnalysisofResults_v2(ptEstimate, concatCCStats, ...
        concatNaiveSigmas, dSigmaCut, segmentDuration, resampleRate1, ...
        window1, window2, figureLegend, [outputFileNamePrefix '_stats.dat'], ...
        1, 1/(1+3/35), outputFileNamePrefix, displayResults, doBadGPSTimes); 
  fprintf('** Done StatisticalAnalysisofResults_v2\n');

  % Figure 8
  [tFFT, omega_t] = FFTofPtEstIntegrand(combinedPtEstInt, resampleRate1, ...
  outputFileNamePrefix, displayResults, figureNumber, ...
  figureLegend);
  fprintf('** Done FFTofPtEstIntegrand\n');


  fprintf('** Done postProcessScriptFull\n');

return; % end postProcessScriptFull
