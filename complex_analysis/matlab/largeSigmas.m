function badGPSTimes = largeSigmas(filename,cutoff);
%
%  function badGPSTimes = largeSigmas(filename,cutoff);
%
%  Reads in theoretical sigmas from a file,
%  and produces a list of "bad" GPS times for which sigmas
%  are too large. 
%  
%  Input:
%
%    filename = name of file with naive and standard theortical sigmas
% 
%    cutoff = minimum allowed value of sigma
%
%    badGPSTimes = column vector with GPS times to veto
%
%  Routine written by V. Mandic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  cutoff;
catch
  cutoff = Inf;
end;

% default return values
badGPSTimes = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in formatted data from a file, extracting relevant input data

[gpsTimes, naiveSigmas, sigmas] = ...
    textread(filename, '%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes) return; end

cc = sigmas > cutoff;
badGPSTimes = gpsTimes(cc);

return
