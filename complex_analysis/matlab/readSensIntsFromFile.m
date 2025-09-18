function [sensInts, flow, deltaF, gpsTimes] = readSensIntsFromFile(filename)
%
%  readSensIntsFromFile -- read sensitivity integrands from file
%
%  [sensInts, flow, deltaF, gpsTimes] = readSensIntsFromFile(filename)
%  returns a 2-d array of sensitivity integrands read in from a file 
%  (integrands corresponding to different times are in different rows).  
%  Also returned are the GPS times and initial frequency and frequency 
%  spacing corresponding to the sensitivity integrands.  The sens ints 
%  are sorted according to GPS time.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: readSensIntsFromFile.m,v 1.1 2005-04-07 14:01:00 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
sensInts = [];
flow = [];
deltaF = [];
gpsTimes = [];

% read in data from file
[gpsTimes, ccSigmas, freqs, sensInts] = ...
    textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes) return; end
                                                                                
[gpsTimes, ind] = unique(gpsTimes);
freqs = unique(freqs);
                                                                                
numFreqs = length(freqs);
flow = freqs(1);
fhigh = freqs(end);
deltaF = freqs(2)-freqs(1);
                                                                                
% create matrix with rows labeled by gps times and columns by frequencies
sensInts = reshape(sensInts, numFreqs, length(gpsTimes));
sensInts = transpose(sensInts);

% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
sensInts = sensInts(ind,:);

return

