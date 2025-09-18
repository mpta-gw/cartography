function [combinedPtEstInt, flow, deltaF, errorBar] = readCombinedPtEstIntFromFile(filename)
%
%  readCombinedPtEstIntFromFile -- read sensitivity integrands from file
%
%  [combinedPtEstInt, flow, deltaF, errorBar] = readCombinedPtEstIntFromFile(filename)
%  returns a combined sensitivity integrand read in from a file
%  Also returned are the initial frequency and frequency spacing 
%  corresponding to the sensitivity integrand, and the error bar 
%  associated with the sensitivity integrand.
%
%  Routine written by Emma L. Robinson
%  Contact elr@star.sr.bham.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default return values
combinedPtEstInt = [];
flow = [];
deltaF = [];
errorBar = [];

% read in data from file
[errorBar, freqs, combinedPtEstIntReal , combinedPtEstIntImag] = ...
    textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(errorBar) return; end
combinedPtEstInt = combinedPtEstIntReal + (1i * combinedPtEstIntImag);

errorBar = errorBar(1);

numFreqs = length(freqs);
flow = freqs(1);
fhigh = freqs(end);
deltaF = freqs(2)-freqs(1);

% create matrix with columns labeled by frequencies
combinedPtEstInt = transpose(combinedPtEstInt);

return


