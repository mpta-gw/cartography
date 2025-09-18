function [ccSpectra, flow, deltaF, gpsTimes] = readCCspectraFromFile(filename)
%
%  readCCspectraFromFile -- read CC spectra from file
%
%  [ccSpectra, flow, deltaF, gpsTimes] = readCCspectraFromFile(filename)
%  returns a 2-d array of CC spectra read in from a file (spectra 
%  corresponding to different times are in different rows).  Also returned
%  are the GPS times and initial frequency and frequency spacing 
%  corresponding to the spectra.  The spectra are sorted according to GPS 
%  time.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: readCCspectraFromFile.m,v 1.1 2005-04-07 14:01:00 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
ccSpectra = [];
flow = [];
deltaF = [];
gpsTimes = [];
 
% read in data from file
[gpsTimes, ccSigmas, freqs, spectra_real, spectra_imag] = ...
  textread(filename, '%f%f%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes) return; end
                                                                                
[gpsTimes, ind] = unique(gpsTimes);
freqs = unique(freqs);
                       
numFreqs = length(freqs);
flow = freqs(1);
fhigh = freqs(end);
deltaF = freqs(2)-freqs(1);

% create matrix with rows labeled by gps times and columns by frequencies
ccSpectra = reshape(spectra_real+1i*spectra_imag, numFreqs, length(gpsTimes));
ccSpectra = transpose(ccSpectra);

% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
ccSpectra = ccSpectra(ind,:);

return

