function prepareCondorRM(paramsFileMat,...
  filePrefix,extraTag,fileSuffix,numJobs,numJobsSuper,...
  segmentDuration, badGPSTimes, ...
  doOverlap,resampleRate1,resampleRate2, hannDuration1, hannDuration2);
%
%  prepareCondorRM --- Used for Radiometer post-processing
%                    put all following parameters into
%                    the .mat paramsFileMatfile.
%                    This makes it easy for each Condor job to load
%                    all relevant parameters
%
%  Routine written by Stefan Ballmer
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try hannDuration1; catch hannDuration1=segmentDuration; end;
try hannDuration2; catch hannDuration2=segmentDuration; end;

if doOverlap
  hannDuration1=segmentDuration;
  hannDuration2=segmentDuration;
end;
numPoints1    = segmentDuration*resampleRate1; 
window1       = tukeywin(numPoints1, hannDuration1/segmentDuration);

numPoints2    = segmentDuration*resampleRate2;
window2       = tukeywin(numPoints2, hannDuration2/segmentDuration);

save(paramsFileMat,'filePrefix','extraTag','fileSuffix','numJobs','numJobsSuper',...
'segmentDuration','badGPSTimes','doOverlap','window1','window2');
return;
