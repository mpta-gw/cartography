function []=createPostProcessingPlanRM(paramsFile,jobsFile,postprocFile)
%
%  createPostProcessingPlanRM --- routine collecting all data required for
%                                 condor post-processing of radiometer data
%  
%  excludeJobs feature not implemented yet.
%
%  Routine written by Stefan Ballmer
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%% read in name/value pairs for post-processing file
[names,values] = ...
  textread(postprocFile, '%s %s\n', -1, 'commentstyle', 'matlab');
%% check that number of names and values are equal
if length(names)~=length(values)
  error('invalid post-processing file');
end
%% loop over parameter names, assigning values to structure
for ii=1:length(names)
  switch names{ii}
    case 'jobsPerSuper'
      jobsPerSuper = str2num(values{ii});
    case 'vetothresh'
      vetothresh = str2num(values{ii});
    case 'combineWhat'
      combineWhat = values{ii};
    case 'excludeJobs'
      excludeJobs = transpose(str2num(values{ii}));
    case 'extraTag'
      extraTag = values{ii};
    case 'restrictSegmentList'
      restrictSegmentList = values{ii};
    case 'paramsFileMat'
      paramsFileMat = values{ii};
    otherwise
      %% do nothing 
  end %% switch
end %% loop over parameter names
% otherwise set default values
try jobsPerSuper;  catch error('jobsPerSuper not defined');  end
try vetothresh;    catch error('vetothresh not defined');    end
try combineWhat;   catch combineWhat='ccstatsSky';           end
try excludeJobs;   catch excludeJobs=[];                     end
try extraTag;      catch extraTag='';                        end
try restSeg=load(restrictSegmentList); catch restSeg=[];     end
try paramsFileMat; catch error('paramsFileMat not defined'); end


% read in job start time and duration from a file
[ignore1, startTimes, ignore2, jobDurations] = ...
  textread(jobsFile, '%n %n %n %n', -1, 'commentstyle', 'matlab'); 

% read in params structure from a file
params = readParamsFromFile(paramsFile);
try params.doDirectional;            catch error('Not a radiometer parameter file (doDirectional missing)'); end;
try params.writeNaiveSigmasToFiles;  catch error('writeNaiveSigmasToFiles missing'); end;
if ~params.doDirectional
  error('Not a radiometer parameter file (doDirectional = false)');
end;
if ~params.writeNaiveSigmasToFiles
  error('Need writeNaiveSigmasToFiles=true for performing veto');
end;
if params.useSkyPatternFile
  params.SkyPattern=load(params.SkyPatternFile);
else
  params.SkyPattern=SkyPattern(params.SkyPatternRightAscensionNumPoints,params.SkyPatternDeclinationNumPoints);
end

			 
filePrefix=[params.outputFilePrefix '_' combineWhat '.job'];
fileSuffix='.trial1.mat';
naivefilePrefix=[params.outputFilePrefix '_naivesigmas.job'];
naivefileSuffix='.trial1.dat';
numJobs=length(startTimes);
segmentDuration=params.segmentDuration;
doOverlap=params.doOverlap;
hannDuration1 = params.hannDuration1;
hannDuration2 = params.hannDuration2;
resampleRate1 = params.resampleRate1;
resampleRate2 = params.resampleRate2;
displayresults = 1;
Nra=params.SkyPatternRightAscensionNumPoints;
Ndecl=params.SkyPatternDeclinationNumPoints;
skyPattern=params.SkyPattern;

for kk=1:numJobs
  filename(kk)={[filePrefix,num2str(kk),fileSuffix]};
  naivname(kk)={[naivefilePrefix,num2str(kk),naivefileSuffix]};
end;



if doOverlap
  hannDuration1=segmentDuration;
  hannDuration2=segmentDuration;
end;
numPoints1    = segmentDuration*resampleRate1; 
window1       = tukeywin(numPoints1, hannDuration1/segmentDuration);

numPoints2    = segmentDuration*resampleRate2;
window2       = tukeywin(numPoints2, hannDuration2/segmentDuration);

badGPSTimes=[];
totSeg=0;
for kk=1:numJobs
  fprintf('Job %d\n',kk);
  fileGood=1;
  try
    p=load(naivname{kk});
  catch p=[]; end;
  lenp=length(p);
  if lenp==0
    fileGood=0;
  end;
  if fileGood
    % the ~( ..< ..) also vetoes of possible NaN's p
    bad= ~(abs(log(p(:,2)./p(:,3))) < log(vetothresh));
    if length(restSeg)>0
      bad=or(bad,~isCompletelyContained(p(:,1),segmentDuration,restSeg));
    end
    badGPSTimes=[badGPSTimes;p(bad,1)];
    totSeg=totSeg+lenp;
  else
    fprintf('Could not load file %s\n',naivname{kk});
  end;
end;
fprintf('Vetoing %d out of %d %d second segments (%f %%)\n',length(badGPSTimes),totSeg,segmentDuration,100.0*length(badGPSTimes)/totSeg);

totJobs=numJobs;
numJobs=1:numJobs;
numJobsSuper=floor((numJobs-1)/jobsPerSuper)+1;
totJobsSuper=numJobsSuper(end);
fprintf('Total jobs / super-jobs: %d / %d\n',totJobs,totJobsSuper);


prepareCondorRM(paramsFileMat,...
  filePrefix,extraTag,fileSuffix,numJobs,numJobsSuper,...
  segmentDuration, badGPSTimes, ...
  doOverlap,resampleRate1,resampleRate2, hannDuration1, hannDuration2);


toc;
