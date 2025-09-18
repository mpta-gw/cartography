function [ccStats, ccSigmas, gpsTimes] ...
      = readCCstatsCCsigmasFromFile(filename, complexflag)
%
%  readCCstatsCCsigmasFromFile -- read CC stats, sigmas from file
%
%  [ccStats, ccSigmas, gpsTimes]
%     = readCCstatsCCsigmasFromFile(filename, complexflag)
%  returns column vectors containing CC statistic values, theoretical 
%  sigmas, and corresponding GPS start times read in from a file.  The 
%  values are sorted according to GPS time.
%
%  If complexflag is false (or omitted), the GPS times, statistic values
%  and theoretical sigmas are assumed to be single columns of real numbers
%
%  If complexflag is true, the statistic values are assumed to be complex
%  with the real and imaginary parts in different columns
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk or john.whelan@ligo.org
%
%  $Id: readCCstatsCCsigmasFromFile.m,v 1.3 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check complex flag
try
  complexflag;
catch
  complexflag = false;
end;

% default return values
ccStats = [];
ccSigmas = [];
gpsTimes = [];

% read in data from file
if complexflag
  [gpsTimes, ccStatsReal, ccStatsImag, ccSigmas] = ...
      textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
  ccStats =  ccStatsReal + 1i*ccStatsImag;
else
  [gpsTimes, ccStats, ccSigmas] = ...
      textread(filename, '%f%f%f\n', -1, 'commentstyle', 'matlab');
end;
if isempty(gpsTimes) return; end
                                                                                
% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
ccStats = ccStats(ind);
ccSigmas = ccSigmas(ind);

return

