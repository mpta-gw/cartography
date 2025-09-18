function ...
  [sky, numSegments] = ...
    combineResultsFromSingleJobRM(jobCell, ...
                      segmentDuration, badGPStimes, ...
                      doOverlap, window1, window2)
%function ...
%  [sky, numSegments] = ...
%    combineResultsFromSingleJobRM(jobCell, ...
%                      segmentDuration, badGPStimes, ...
%                      doOverlap, window1, window2)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands.
%
%  Input:
%
%    jobCell - cell array (1 entry per segment)
%              each entry is a struct
%               - time: GPS time
%               - data: Nx2 array, 1st column = ccStat, 2nd column = sigma
%                       rows correspond to different points in the sky
%              can also be a filename
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
%  $Id: combineResultsFromSingleJobRM.m,v 1.2 2005-09-28 20:38:10 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(iscell(jobCell))
  try
    p=load(jobCell);
    jobCell=p.Sky;
  catch
    try
      fname=[jobCell(1:end-4),'.dat'];
      [jobCell,x]=ccDataLoad(fname);
    catch
      fprintf('Could not load file %s\n',jobCell);
      sky=[];
      numSegments=0;
      return;
    end;
  end;
end

numSegmentsTotal = length(jobCell);
numSegments = 0;
if numSegmentsTotal ==0
  sky=[];
  return;
end

if doOverlap
  stepsize=segmentDuration/2;
else
  stepsize=segmentDuration;
end

ii=1;
first=1;
last=0;
expectedTime=jobCell{1}.time;
skycount=1;

  while ii <= numSegmentsTotal
    bad=ismember(jobCell{ii}.time,badGPStimes);
    time=jobCell{ii}.time;
    if or(bad,time~=expectedTime)
      if first<=last  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate
	  %fprintf('%d\t%d\n',jobCell{first}.time,jobCell{last}.time);
          [map{skycount}.data, numSegmentsPart] = ...
            combineResultsRM(jobCell, ...
                             segmentDuration, first, last, ...
                             doOverlap, window1, window2);
	  sky{skycount}.data=jobCell{first}.time;
	  skycount=skycount+1;
	  if numSegmentsPart ~= last-first+1
	    error('combineResultsRM failed');
	    return;
	  end
	  numSegments=numSegments + last-first+1;
	first=ii; last=ii;
	if bad
	  first=first+1;
	end
      else  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if bad
	  first=first+1;
	end
        last=last+1;
      end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      last=last+1;
    end
    ii=ii+1;
    expectedTime=time+stepsize;
  end
      if first<=last
        % calculate
	  %fprintf('%d\t%d\n',jobCell{first}.time,jobCell{last}.time);
          [map{skycount}.data, numSegmentsPart] = ...
            combineResultsRM(jobCell, ...
                             segmentDuration, first, last, ...
                             doOverlap, window1, window2);
	  sky{skycount}.data=jobCell{first}.time;
	  skycount=skycount+1;
	  if numSegmentsPart ~= last-first+1
	    error('combineResultsRM failed');
	    return;
	  end
	  numSegments=numSegments + last-first+1;
      end
  if numSegments == 0
    sky=[];
  elseif skycount>2
    [sky, numSegmentsPart] = ...
            combineResultsRM(map, ...
                             1, 1, skycount-1, ...
                             false, window1, window2);
  else
    sky=map{1}.data;
  end;
return

