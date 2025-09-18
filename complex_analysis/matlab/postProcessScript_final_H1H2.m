function postProcessScript_final_H1H2(postParamsfile,varargin)
%    postProcessScript_final_H1H2(postParamfile,jobNumber)
%    postProcessScript_final_H1H2(postParamfile,startGPS,endGPS)
% This program postprocess the files obtained from running stochastic.m .
% The input for this code has either two or three arguments. 
% The first argument is the name of a file which contains a few necessary 
% parameters as name/value pairs (e.g., postParamfile.txt). It is similar 
% to the parameter file used for running stochastic.m, but has less number 
% of parameters (and a few parameters different from the earlier one).
%    we can either provide a job number (like 123) or a two component vector 
% specifying initial job number and final job number (like [123 456]) or start 
% and end GPS Times (like 815413494,815565790). If only a single job number 
% provided and 'runSingleJob' field in the paramfile kept 'true' it will 
% only run that job otherwise it will run all the jobs from 1 to that 
% job number. The routine can also postprocess the results with different 
% frequency mask, if requested (Note that it can notch only additional 
% frequencies and cannot be used for a entirely new frequency notching). 
% 
% This can also be compiled and run (but two component vector input for the 
% job number doesn't work ).
%
% For any suggestions and comments
% Contact shivara@physics.umn.edu 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin==2)
   jobNumber = strassign(varargin{1});% required for compiled code
   runFlag = 0; % 1 for GPS, 0 for Job nos
elseif(nargin==3) 
   StarT = strassign(varargin{1});
   EnD = strassign(varargin{2});
   runFlag = 1; % 1 for GPS, 0 for Job nos
else
   error('I am sorry, the no. of arguments should be 2 or 3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in necessary name/value pairs using ctextread and assigning
% them to varaibles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[names,values] = ...
  ctextread(postParamsfile, '%s %s\n', -1, 'commentstyle', 'matlab');

% check that the number of names and values are equal
if length(names)~=length(values)
  error('invalid parameter file');
end

% loop over parameter names, assigning values to structure
for ii=1:length(names)
 switch names{ii}
   case 'stoc_files'
     stoc_files = values{ii};
   case 'jobFile'
     jobsFile = values{ii};
   case 'runSingleJob'
     runSingleJob = str2num(values{ii});
   case 'CombinedFiles'
     CombinedFiles = str2num(values{ii});
   case 'badGPSTimes'
     badGPSTimes = values{ii};
   case 'flow'
     flow = str2num(values{ii});
   case 'fhigh'
     fhigh = str2num(values{ii});
   case 'figureLegend'
     figureLegend = values{ii};
   case 'outputFileNamePrefix'
     outputFileNamePrefix = values{ii};
   case 'doOverlap'
     doOverlap = str2num(values{ii});
   case 'segmentDuration'
     segmentDuration = str2num(values{ii});
   case 'resampleRate'
     resampleRate = str2num(values{ii});
   case 'deltaF'
     deltaF = str2num(values{ii});
   case 'displayresults'
     displayresults = str2num(values{ii});
   case 'addlFreqNotch'
     addlFreqNotch = str2num(values{ii});
   case 'addlFreqRemove'
     addlFreqRemove = transpose(str2num(values{ii}));
   case 'addlBinRemove'
     addlBinRemove = transpose(str2num(values{ii}));
   otherwise
     % do nothing
 end % switch   
end % for loop over ii
  
if(~CombinedFiles)
 commonFilePrefix=[stoc_files];
else
 commonFilePrefix=[stoc_files '_combined'];
end 

if(~exist('complexflag','var'))
  complexflag = 0;
end

% checking for badGPSTimes file 
if(exist('badGPSTimes')~=0 && CombinedFiles==0)
 badGPSTimes = load(badGPSTimes);
else
 badGPSTimes=[];
end

% New output file name and legend which includes start and end GPS or job No
if(nargin == 3) 
 outputFileNamePrefix = [outputFileNamePrefix '_' num2str(StarT) '_' num2str(EnD)];
 figureLegend = [figureLegend ' ' num2str(StarT) ' '  num2str(EnD)];
elseif(nargin == 2 && runSingleJob == 1 && length(jobNumber)==1)
 outputFileNamePrefix = [outputFileNamePrefix '_job' num2str(jobNumber)];
 figureLegend = [figureLegend ' jobNo ' num2str(jobNumber)];
end

% constructing additional freq. mask, if requested.
if addlFreqNotch
 modifyFilterC = constructFreqMask(flow, fhigh, deltaF, addlFreqRemove, addlBinRemove, addlFreqNotch);
 modifyFilter = transpose(modifyFilterC);
else
 modifyFilter = 1;
end

% Check whether combinedFiles and addlFreqNotch true at same time; 
% addlFreqNotch works only for individual ccstat, ccspectra  and sensint files
if (addlFreqNotch & CombinedFiles)
  error('addlFreqNotch and CombinedFiles both are true');
end


% Deciding the start and end job nos, in case the input was job numbers
if ~runFlag
 if(length(jobNumber)==2)
  StarT = jobNumber(1);
  EnD = jobNumber(2);
 elseif(length(jobNumber) == 1 && runSingleJob == 1)
  StarT = jobNumber;
  EnD = jobNumber;
 else
  StarT = 1;
  EnD = jobNumber;
 end
end

% some variables common to all the analyses
window1=hann(segmentDuration*resampleRate);
window2=hann(segmentDuration*resampleRate);
fileSuffix='.trial1.dat';
figureNumber = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combining the results for the requested Job(s) or period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ptEstimate, errorBar, combinedPtEstInt, combinedSensInt, ...
          numSegmentsTotal] = ...
  combineResultsFromMultipleJobs_H1H2(commonFilePrefix, fileSuffix, ...
                                 StarT, EnD, segmentDuration, badGPSTimes, ...
                                 addlFreqNotch, modifyFilter, ...
                                 doOverlap, window1, window2, ...
                                 outputFileNamePrefix, displayresults, ...
                                 figureNumber, figureLegend, ...
                                 runFlag, jobsFile, complexflag);

if addlFreqNotch
  outputFileNamePrefix = [outputFileNamePrefix '_renormalized'];
end

if combinedPtEstInt.data == 0
  return;
end

[tFFT, omega_t] = FFTofPtEstIntegrand(combinedPtEstInt, resampleRate, outputFileNamePrefix, displayresults, 21, figureLegend);

return
