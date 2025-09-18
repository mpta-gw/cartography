function ...
  [sky, numSegments] = ...
    combineResultsFromMultipleJobsRM(filePrefix,fileSuffix,jobNums, ...
                      segmentDuration, badGPStimes, ...
                      doOverlap, window1, window2)
%function ...
%  [sky, numSegments] = ...
%    combineResultsFromMultipleJobsRM(filePrefix,fileSuffix,jobNums, ...
%                      segmentDuration, badGPStimes, ...
%                      doOverlap, window1, window2)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands.
%
%  Input:
%
%    filePrefix  - common filename prefix
%    fileSuffix  - common filename suffix
%    jobNums     - Vector with job numbers
%    segmentDuration = length of analysis segment in sec (typically 60)
%    badGPStimes - vector of GPS segment start times to be ignored
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%
%  Output:
%
%    sky      = Nx2 array with optimal point estimate of Omega0 (column 1) and 
%              theoretical error bar for the point estimate (column 2)
%              rows correspond to different points in the sky
%    numSegments = number of data segments combined
%
%  Routine written by Stefan Ballmer
%  Contact sballmer@ligo.mit.edu
%
%  $Id: combineResultsFromMultipleJobsRM.m,v 1.2 2005-08-11 15:53:59 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSegments=0;
skycount=1;
map.time=0;

for kk=jobNums
  filename=[filePrefix,num2str(kk),fileSuffix];
  [map.data, numSegmentsPart] = ...
    combineResultsFromSingleJobRM(filename, ...
                      segmentDuration, badGPStimes, ...
                      doOverlap, window1, window2);
  if numSegmentsPart>0
    jobmap{skycount}=map;
    skycount=skycount+1;
    numSegments=numSegments+numSegmentsPart;
  end
end

  if numSegments == 0
    sky=[];
  elseif skycount>2
    [sky, numSegmentsPart] = ...
            combineResultsRM(jobmap, ...
                             1, 1, skycount-1, ...
                             false, window1, window2);
  else
    sky=jobmap{1}.data;
  end;

return

