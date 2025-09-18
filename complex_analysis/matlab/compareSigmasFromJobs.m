function [gpsTimes, sigmas, naiveSigmas, badGPSTimes, ...
	  goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
      compareSigmasFromJobs(filePrefix, fileSuffix, numJobs, ...
			    minRatio, maxRatio);
%
%  [gpsTimes, sigmas, naiveSigmas, badGPSTimes, ...
% 	  goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
%       compareSigmasFromJobs(filePrefix, fileSuffix, numJobs, ...
% 			    minRatio, maxRatio);
%
%  Reads in "naive" and actual theoretical sigmas from a collection of
%  jobs, by calling 'compareSigmasFromFile' repeatedly,
%  and produces a list of "bad" GPS times for which the ratio of
%  the two sigmas is above or below specified threshholds
%  
%  Input:
%
%    filePrefix, fileSuffix: prefix, suffix of filenames containing
%    the "naive" and actual theoretical sigmas
%    (note: the filename is assumed to be of the form
%    filePrefixNfileSuffix, where N runs from 1 to numJobs
%
%    numJobs = total number of jobs 
%
%    minRatio = minimum "acceptable" ratio of (theor sigma)/(naive theor sigma)
%
%    maxRatio = maximum "acceptable" ratio of (theor sigma)/(naive theor sigma)
%
%  Output:
%
%    gpsTimes = column vector with GPS times of all jobs
%
%    sigmas = column vector with theoretical sigma values
%
%    naiveSigmas = column vector with "naive" theoretical sigma values
%
%    badGPSTimes = column vector with GPS times to veto
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: compareSigmasFromJobs.m,v 1.5 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize arrays
gpsTimes = [];
sigmas = [];
naiveSigmas = [];
badGPSTimes = [];
goodGPSTimes = [];
goodSigmas = [];
goodNaiveSigmas = [];

% loop over number of jobs
for k = 1:numJobs
 
%  fprintf('Analysing job %d\n', k);

  filename = [filePrefix num2str(k) fileSuffix];

  % check for missing job files
  if exist(filename)
     [myGPSTimes, mySigmas, myNaiveSigmas, myBadGPSTimes, ...
      myGoodGPSTimes, myGoodSigmas, myGoodNaiveSigmas] = ...
     compareSigmasFromFile(filename, minRatio, maxRatio);

     if size(myGPSTimes)==0
        fprintf('Job %d does not contain any cc spectra\n', k);
     else
        gpsTimes = [gpsTimes; myGPSTimes];
        sigmas = [sigmas; mySigmas];
        naiveSigmas = [naiveSigmas; myNaiveSigmas];
        badGPSTimes = [badGPSTimes; myBadGPSTimes];
        goodGPSTimes = [goodGPSTimes; myGoodGPSTimes];
        goodSigmas = [goodSigmas; myGoodSigmas];
        goodNaiveSigmas = [goodNaiveSigmas; myGoodNaiveSigmas];
     end;

  else
    fprintf('Missing job %d (filename %s)\n', k, filename)
  
  end

end

return
