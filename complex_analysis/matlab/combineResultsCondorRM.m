function combineResultsCondorRM(paramsFileMat,superJobNum)
%
%  combineResultsCondorRM --- Runs one super-job combining the output
%                             radiometer maps from several jobs
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

superJobNum=strassign(superJobNum);
p=load(paramsFileMat);
ind=find(p.numJobsSuper==superJobNum);
Sky{1}.time=0;
[Sky{1}.data, numSegments] = combineResultsFromMultipleJobsRM(p.filePrefix,p.fileSuffix,p.numJobs(ind),...
                                                      p.segmentDuration, p.badGPSTimes, ...
						      p.doOverlap, p.window1, p.window2);
outfile=[p.filePrefix,p.extraTag,'Super',num2str(superJobNum),p.fileSuffix];
save(outfile,'Sky');
return
