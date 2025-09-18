% script to produce a single file containing a list of bad GPS times
%
% concatenates results from abs sigma cut, airplane veto, acoustic veto,
% seismic veto, and S3 stochastic hardware injection times.
%
% $Id: badGPSTimes.m,v 1.1 2005-10-19 17:47:30 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ddmmyyyyhhmmss  = datestr(now);

if 1% H2-L1 hann
sigmasfile='../output/ccsigmas_H2L1_hann_sorted.dat';
abssigmaFile = '../output/H2L1_hann_vetotimes_abssigmacut.dat';
largesigmaFile = '../output/H2L1_hann_vetotimes_largesigmacut.dat';
ccStatsFile  = '../output/ccstats_H2L1_hann_sorted.dat';
acousticFile = '';
airplaneFile = '';
%seismicFile  = '';
seismicFile  = '../output/H2L1seismic.dat';
H1hardwareFile = '';
H2hardwareFile = '';
L1hardwareFile = '';
outputfilename = '../output/badGPSTimes_H2L1_hann.dat';
minRatio = 0.8;
maxRatio = 1.2;
cutoff = 1;
end

if 0% H1-L1 hann
sigmasfile='../output/ccsigmas_H1L1_hann_sorted.dat';
abssigmaFile = '../output/H1L1_hann_vetotimes_abssigmacut_2.dat';
largesigmaFile = '../output/H1L1_hann_vetotimes_largesigmacut.dat';
ccStatsFile  = '../output/ccstats_H1L1_hann_sorted.dat';
acousticFile = '';
airplaneFile = '';
%seismicFile  = '';
seismicFile  = '../output/H1L1seismic.dat';
H1hardwareFile = '';
H2hardwareFile = '';
L1hardwareFile = '';
outputfilename = '../output/badGPSTimes_H1L1_hann_3.dat';
minRatio = 0.8;
maxRatio = 1.2;
cutoff = 1;
end

if 0% H1-L1 hann
sigmasfile='../output/ccsigmas_H1L1_hann_sorted.dat';
abssigmaFile = '../output/H1L1_hann_vetotimes_abssigmacut_2.dat';
ccStatsFile  = '../output/ccstats_H1L1_hann_sorted.dat';
acousticFile = '';
airplaneFile = '';
%seismicFile  = '';
seismicFile  = '../output/H1L1seismic.dat';
H1hardwareFile = '';
H2hardwareFile = '';
L1hardwareFile = '';
outputfilename = '../output/badGPSTimes_H1L1_hann_2.dat';
minRatio = 0.8;
maxRatio = 1.2;
largesigmaFile = '';
cutoff = 0;
end

if 0% H1-L1 hann
sigmasfile='../output/ccsigmas_H1L1_hann_sorted.dat';
abssigmaFile = '../output/H1L1_hann_vetotimes_abssigmacut.dat';
ccStatsFile  = '../output/ccstats_H1L1_hann_sorted.dat';
acousticFile = '';
airplaneFile = '';
seismicFile  = '';
H1hardwareFile = '';
H2hardwareFile = '';
L1hardwareFile = '';
outputfilename = '../output/badGPSTimes_H1L1_hann.dat';
minRatio = 0.8;
maxRatio = 1.2;
largesigmaFile = '';
cutoff = 0;
end

if 0% H1-H2 hann
sigmasfile='../output/ccsigmas_H1H2_hann_sorted.dat';
abssigmaFile = '../output/H1H2_hann_vetotimes_abssigmacut.dat';
ccStatsFile  = '../output/ccstats_H1H2_hann_sorted.dat';
acousticFile = '';
airplaneFile = '';
seismicFile  = '';
H1hardwareFile = '';
H2hardwareFile = '';
L1hardwareFile = '';
outputfilename = '../output/badGPSTimes_H1H2_hann.dat';
minRatio = 0.8;
maxRatio = 1.2;
largesigmaFile = '';
cutoff = 0;
end


[gpsTimes, sigmas, naiveSigmas, badTimes_abssigma, ...
  goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
      compareSigmasFromFile(sigmasfile, minRatio, maxRatio);
fprintf('Total number of segments = %d\n',length(gpsTimes));
fprintf('Number of good segments  = %d\n',length(goodSigmas));

% write bad times to file
fid = fopen(abssigmaFile, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\n', '%gps sec');
for ii=1:length(badTimes_abssigma)
  fprintf(fid, '%d\n', badTimes_abssigma(ii));
end
fclose(fid);


if ~isempty(largesigmaFile)
  badTimes_largesigma = largeSigmas(sigmasfile,cutoff);
  
  % write bad times to file
  fid = fopen(largesigmaFile, 'w');
  fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
  fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
  fprintf(fid, '%s\n', '%gps sec');
  for ii=1:length(badTimes_largesigma)
    fprintf(fid, '%d\n', badTimes_largesigma(ii));
  end
  fclose(fid);
end

% common variables
segmentDuration=60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise arrays
badTimes_acoustic=[];
badTimes_airplane=[];
badTimes_seismic=[];
badTimes_H1hardware=[];
badTimes_H2hardware=[];
badTimes_L1hardware=[];

% initialise counters
counter1 = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;
counter5 = 0;
counter6 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in original analysed times
data = load(ccStatsFile,'ascii'); 
times_start = data(:,1);
times_stop = data(:,1);

% read in acoustic veto times
if ~isempty(acousticFile)
  data = load(acousticFile);
  acoustic_start = data(:,1);
  acoustic_stop  = data(:,2);
end;

% read in airplane veto times
if ~isempty(airplaneFile)
  data = load(airplaneFile);
  airplane_start = data(:,1);
  airplane_stop  = data(:,2);
end;

% read in seismic veto times
if ~isempty(seismicFile)
  data = load(seismicFile);
  seismic_start = data(:,1);
  seismic_stop  = data(:,2);
end;

% read in H1 hardware injection veto times
if ~isempty(H1hardwareFile)
  data = load(H1hardwareFile);
  H1hardware_start = data(:,1);
  H1hardware_stop  = data(:,2);
end;

% read in H2 hardware injection veto times
if ~isempty(H2hardwareFile)
  data = load(H2hardwareFile);
  H2hardware_start = data(:,1);
  H2hardware_stop  = data(:,2);
end;

% read in L1 hardware injection veto times
if ~isempty(L1hardwareFile)
  data = load(L1hardwareFile);
  L1hardware_start = data(:,1);
  L1hardware_stop  = data(:,2);
end;

% loop over original analysed times
for ii=1:length(times_start)

  % acoustic veto
  if ~isempty(acousticFile)
    for jj=1:length(acoustic_start)
      if (  ( (acoustic_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= acoustic_stop(jj)) ) |... 
            ( (acoustic_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= acoustic_stop(jj) ) ) )
        counter1 = counter1+1;
        badTimes_acoustic(counter1) = times_start(ii);
      end
    end
  end

  % airplane veto
  if ~isempty(airplaneFile)
    for jj=1:length(airplane_start)
      if (  ( (airplane_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= airplane_stop(jj)) ) |... 
            ( (airplane_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= airplane_stop(jj) ) ) )
        counter2 = counter2+1;
        badTimes_airplane(counter2) = times_start(ii);
      end
    end
  end

  % seismic veto
  if ~isempty(seismicFile)
    for jj=1:length(seismic_start)
      if (  ( (seismic_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= seismic_stop(jj)) ) |... 
            ( (seismic_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= seismic_stop(jj) ) ) )
        counter3 = counter3+1;
        badTimes_seismic(counter3) = times_start(ii);
      end
    end
  end

  % H1 hardware veto
  if ~isempty(H1hardwareFile)
    for jj=1:length(H1hardware_start)
      if (  ( (H1hardware_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= H1hardware_stop(jj)) ) |... 
            ( (H1hardware_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= H1hardware_stop(jj) ) ) )
        counter4 = counter4+1;
        badTimes_H1hardware(counter4) = times_start(ii);
      end
    end
  end

  % H2 hardware veto
  if ~isempty(H2hardwareFile)
    for jj=1:length(H2hardware_start)
      if (  ( (H2hardware_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= H2hardware_stop(jj)) ) |... 
            ( (H2hardware_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= H2hardware_stop(jj) ) ) )
        counter5 = counter5+1;
        badTimes_H2hardware(counter5) = times_start(ii);
      end
    end
  end

  % L1 hardware veto
  if ~isempty(L1hardwareFile)
    for jj=1:length(L1hardware_start)
      if (  ( (L1hardware_start(jj) <= times_start(ii)) &...
              (times_start(ii) <= L1hardware_stop(jj)) ) |... 
            ( (L1hardware_start(jj) <= times_stop(ii)) &...
              (times_stop(ii) <= L1hardware_stop(jj) ) ) )
        counter6 = counter6+1;
        badTimes_L1hardware(counter6) = times_start(ii);
      end
    end
  end

end % loop over ii (times_tot)

% read in abssigma veto times 
% (these times are already in the desired form)
if ~isempty(abssigmaFile)
  badTimes_abssigmacut = load(abssigmaFile,'ascii'); 
end

% read in large-sigma veto times 
% (these times are already in the desired form)
if ~isempty(largesigmaFile)
  badTimes_largesigmacut = load(largesigmaFile,'ascii'); 
else
  badTimes_largesigmacut = []; 
end

% add on bad abs sigma cut times
badTimes = [transpose(badTimes_acoustic); transpose(badTimes_airplane); ...
            transpose(badTimes_seismic); transpose(badTimes_H1hardware); ...
            transpose(badTimes_H2hardware); transpose(badTimes_L1hardware); ...
            badTimes_abssigmacut; badTimes_largesigmacut];

badTimes = sort(unique(badTimes));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display fraction lost to various cuts
fprintf('Total number of segments analysed (before cuts) = %d\n', length(times_start));
fprintf('Number of segments lost to acoustic veto = %d, fraction = %g\n', length(badTimes_acoustic), length(badTimes_acoustic)/length(times_start));
fprintf('Number of segments lost to airplane veto = %d, fraction = %g\n', length(badTimes_airplane), length(badTimes_airplane)/length(times_start));
fprintf('Number of segments lost to seismic veto = %d, fraction = %g\n', length(badTimes_seismic), length(badTimes_seismic)/length(times_start));
fprintf('Number of segments lost to H1 hardware injection = %d, fraction = %g\n', length(badTimes_H1hardware), length(badTimes_H1hardware)/length(times_start));
fprintf('Number of segments lost to H2 hardware injection = %d, fraction = %g\n', length(badTimes_H2hardware), length(badTimes_H2hardware)/length(times_start));
fprintf('Number of segments lost to L1 hardware injection = %d, fraction = %g\n', length(badTimes_L1hardware), length(badTimes_L1hardware)/length(times_start));
fprintf('Number of segments lost to abs sigma cut = %d, fraction = %g\n', length(badTimes_abssigmacut), length(badTimes_abssigmacut)/length(times_start));
fprintf('Number of segments lost to large sigma cut = %d, fraction = %g\n', length(badTimes_largesigmacut), length(badTimes_largesigmacut)/length(times_start));
fprintf('Total number of segments lost to all vetoes = %d, fraction = %g\n', length(badTimes), length(badTimes)/length(times_start));

% write bad GPS times to file
fid = fopen(outputfilename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', ddmmyyyyhhmmss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\n', '%gps sec');
for ii=1:length(badTimes)
  fprintf(fid, '%d\n', badTimes(ii));
end
fclose(fid);


