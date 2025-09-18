function [optFilts, flow, deltaF, gpsTimes] = readOptFiltsFromFile(filename)
%
%  readOptFiltsFromFile -- read optimal filters from file
%
%  [optFilts, flow, deltaF, gpsTimes] = readOptFiltsFromFile(filename)
%  returns a 2-d array of optimal filters read in from a file 
%  (integrands corresponding to different times are in different rows).  
%  Also returned are the GPS times and initial frequency and frequency 
%  spacing corresponding to the optimal filters.  The sens ints 
%  are sorted according to GPS time.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: readOptFiltsFromFile.m,v 1.2 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
optFilts = [];
flow = [];
deltaF = [];
gpsTimes = [];

% read in data from file
[gpsTimes, freqs, optFilts_real, optFilts_imag] = ...
    textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes) return; end

[gpsTimes, ind] = unique(gpsTimes);
freqs = unique(freqs);

numFreqs = length(freqs);
flow = freqs(1);
fhigh = freqs(end);
deltaF = freqs(2)-freqs(1);

% create matrix with rows labeled by gps times and columns by frequencies
optFilts = reshape(optFilts_real+1i*optFilts_imag, numFreqs, length(gpsTimes));
optFilts = transpose(optFilts);

% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
optFilts = optFilts(ind,:);

return
