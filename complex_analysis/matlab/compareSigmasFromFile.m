function [gpsTimes, sigmas, naiveSigmas, badGPSTimes, ...
	  goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
      compareSigmasFromFile(filename, minRatio, maxRatio);
%
%  [gpsTimes, sigmas, naiveSigmas, badGPSTimes, ...
%	  goodGPSTimes, goodSigmas, goodNaiveSigmas] = ...
%      compareSigmasFromFile(filename, minRatio, maxRatio);
%
%  Reads in "naive" and actual theoretical sigmas from a file,
%  and produces a list of "bad" GPS times for which the ratio of
%  the two sigmas is above or below specified threshholds
%  
%  Input:
%
%    filename = name of file with naive and standard theortical sigmas
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
%  $Id: compareSigmasFromFile.m,v 1.4 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  minRatio;
catch
  minRatio = 0;
end;

try
  maxRatio;
catch
  maxRatio = Inf;
end;

% default return values
gpsTimes = [];
sigmas = [];
naiveSigmas = [];
badGPSTimes = [];
goodGPSTimes = [];
goodSigmas = [];
goodNaiveSigmas = [];
                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in formatted data from a file, extracting relevant input data

[gpsTimes, naiveSigmas, sigmas] = ...
    textread(filename, '%f%f%f\n', -1, 'commentstyle', 'matlab');
  if isempty(gpsTimes) return; end

ratios = sigmas ./ naiveSigmas;

badInd = find( (ratios > maxRatio) | (ratios < minRatio) );
goodInd =  find( (ratios <= maxRatio) & (ratios >= minRatio) );

badGPSTimes = gpsTimes(badInd);

goodGPSTimes =  gpsTimes(goodInd);

goodSigmas = sigmas(goodInd);

goodNaiveSigmas = naiveSigmas(goodInd);

return
