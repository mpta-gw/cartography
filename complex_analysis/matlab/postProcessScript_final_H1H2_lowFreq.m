function postProcessScript_final_H1H2_lowFreq(postParamsfile,varargin)
%   postProcessScript_final_H1H2(postParamfile,jobNumber,addl_freq_mask)
%   postProcessScript_final_H1H2(postParamfile,startGPS,endGPS,addl_freq_mask)
% This program postprocess the files obtained from running stochastic.m .
% The input for this code has either three or four arguments. 
% The first argument is the name of a file which contains a few necessary 
% parameters as name/value pairs (e.g., postParamfile.txt). It is similar 
% to the parameter file used for running stochastic.m, but has less number 
% of parameters (and a few parameters different from the earlier one).
%    The parameter jobNumber can be either a single number (like 123) or a two 
% component vector specifying initial job number and final job number 
% (like [123 456]). If only a single job number is provided and 'runSingleJob' 
% field in the paramfile is kept 'true' then it will run only that particular 
% job or otherwise (if runSingleJob='false' or it is omitted from the 
% postParamsfile) it will run all the jobs from 1 to that job number.
%    Instead of job number, if we want, we can also use startGPS and endGPS 
% which will be of the standard form (like 815413494). 
%
%     The routine can also postprocess the results with different 
% frequency mask, if requested (note that it can notch only additional 
% frequencies that were not notched in the initial stage and cannot be used for
% a entirely new frequency notching). The additional frequency mask should be
% a two column matrix with the first column corresponding to frequencies to be
% notched while the second column representing the no of bins to be notched.
% It should be noted that this argument mondatory, so if we don't want to 
% notch any new frequencies then we should pass an empty matrix [].  
% 
% This code can be compiled and run (but two component vector input for the 
% job number doesn't work ).
%
% For any suggestions and comments
% Contact shivara@physics.umn.edu 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin==3)
   jobNumber = strassign(varargin{1});% required for compiled code
   runFlag = 0; % 1 for GPS, 0 for Job nos
   remove_freq = varargin{2};
elseif(nargin==4) 
   StarT = strassign(varargin{1});
   EnD = strassign(varargin{2});
   runFlag = 1; % 1 for GPS, 0 for Job nos
   remove_freq = varargin{3};
else
   error('Sorry, the no. of arguments should be 3 or 4');
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
   case 'fileSuffix'
     fileSuffix = values{ii};
%   case 'addlFreqRemove'
%     addlFreqRemove = transpose(str2num(values{ii}))
%   case 'addlBinRemove'
%     addlBinRemove = transpose(str2num(values{ii}))
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
if(nargin == 4) 
 outputFileNamePrefix = [outputFileNamePrefix '_' num2str(StarT) '_' num2str(EnD)];
 figureLegend = [figureLegend ' ' num2str(StarT) ' '  num2str(EnD)];
elseif(nargin == 3 && runSingleJob == 1 && length(jobNumber)==1)
 outputFileNamePrefix = [outputFileNamePrefix '_job' num2str(jobNumber)];
 figureLegend = [figureLegend ' jobNo ' num2str(jobNumber)];
end

% constructing additional freq. mask, if requested.
if(addlFreqNotch & ~isempty(remove_freq))
 addlFreqRemove = remove_freq(:,1);
 addlBinRemove = remove_freq(:,2);
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
figureNumber = 10; % just a number

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
