function [y, gpsTimes] = ignoreBadGPSTimes(x, totalGPSTimes, badGPSTimes)
%
%  ignoreBadGPSTimes  -- ignores data corresponding to bad GPS times
%
%  [y, gpsTimes] = ignoreBadGPSTimes(x, totalGPSTimes, badGPSTimes)
%  returns an array of good GPS times and corresponding good data values
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: ignoreBadGPSTimes.m,v 1.1 2005-04-14 12:45:19 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gpsTimes, ind] = setdiff(totalGPSTimes, badGPSTimes);
y = x(ind,:);

return

