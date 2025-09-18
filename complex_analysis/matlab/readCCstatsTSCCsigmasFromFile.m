function [ccStats, ccSigmas, tau, gpsTimes] = readCCstatsTSCCsigmasFromFile(filename)
%
%  readCCstatsTSCCsigmasFromFile -- read CC spectra from file
%
%  [ccStats, ccSigmas, tau, gpsTimes] = readCCstatsTSCCsigmasFromFile(filename)
%  returns a 2-d array of CC stat time series read in from a file.  Also returned
%  are the GPS times and time shift vector tau. The correlation time series
%  are sorted according to GPS  time.
%
%  Routine written by Joseph D. Romano / Stefan Ballmer.
%  Contact Joseph.Romano@astro.cf.ac.uk / sballmer@ligo.mit.edu
%
%  $Id: readCCstatsTSCCsigmasFromFile.m,v 1.1 2005-05-19 23:01:46 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
ccStats = [];
ccSigmas = [];
tau = [];
gpsTimes = [];
 
% read in data from file
[gpsTimes, ccSigmas, tau, ccStats] = ...
  ctextread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes) return; end
                                                                                
[gpsTimes, ind] = unique(gpsTimes);
tau = unique(tau);
                       
numTimes = length(tau);
%flow = freqs(1);
%fhigh = freqs(end);
%deltaF = freqs(2)-freqs(1);

% create matrix with rows labeled by gps times and columns by frequencies
ccStats = reshape(ccStats, numTimes, length(gpsTimes));
ccStats = transpose(ccStats);

% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
ccStats = ccStats(ind,:);

return

