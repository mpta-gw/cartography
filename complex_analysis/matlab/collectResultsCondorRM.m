function collectResultsCondorRM(paramsFileMat)
%
%  collectResultsCondorRM --- Collects the radiometer data from several
%                             super-jobs.
%                             All relevant parameters are passed through
%                             the .mat file paramsFileMat, created by
%                             prepareCondor.
%
%    paramsfile is a mat file containing the following variables:
%       filePrefix,fileSuffix,numJobs,numJobsSuper,...
%       segmentDuration, badGPSTimes, ...
%       doOverlap, window1, window2
%
%  Routine written by Stefan Ballmer
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=load(paramsFileMat);
filePrefixS=[p.filePrefix,p.extraTag,'Super'];

[map, numSegments] = combineResultsFromMultipleJobsRM(filePrefixS,p.fileSuffix, unique(p.numJobsSuper),...
                                                      1, [], ...
						      false, p.window1, p.window2);
outfile=[p.filePrefix,p.extraTag,'FINALMAP',p.fileSuffix];
save(outfile,'map');
return
