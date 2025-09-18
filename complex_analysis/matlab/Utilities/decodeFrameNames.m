function [gpsTimes,frameDurs]=decodeFrameNames(frameFiles)
%
%  decodeFrameNames --- extract GPS times and frame durations from names
%
%  [gpsTimes,frameDurs]=decodeFrameNames('frameFiles') parses frames
%  with the naming convention <site>-<type>-GPStime-duration.gwf as
%  specified in http://www.ligo.caltech.edu/docs/T/T010150-00.pdf
%
%  The input is
%
%   frameFiles  -  a cell array of frame file names.
%
%  The outputs are
%
%   gpsTimes  - vector of gps start times
%   frameDurs  - vector of frame file durations
%
%  Routine cannibalized from Keith Thorne's dir2framelist.m in the
%  matapps channel package by John Whelan
%  Contact john.whelan@ligo.org
%
% $Id: decodeFrameNames.m,v 1.2 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nFiles = length(frameFiles);
    gpsList = zeros(nFiles,1);
    durList = zeros(nFiles,1);
    for iFile = 1:nFiles
        thisFile = frameFiles{iFile,1};
        lastChar = length(thisFile);
        dashPos = strfind(thisFile,'-');
        numDash = numel(dashPos);
        if(numDash > 1)
            gpsBeg = dashPos(numDash-1) + 1;
            gpsEnd = dashPos(numDash) - 1;
            gpsStart = str2num(thisFile(gpsBeg:gpsEnd));
            durBeg = dashPos(numDash)+1;
            durEnd = lastChar - 4;
            durVal = str2num(thisFile(durBeg:durEnd));
        else
            msgId = 'decodeFrameNames:badName';
            warning(msgId,'%s: file %s not to frame-file name convention',...
                msgid,fileList(iFile).name);
            gpsStart = 0;
            gpsEnd = 0;
        end
        gpsList(iFile) = gpsStart;
        durList(iFile) = durVal;
    end
    [gpsTimes,ndx] = unique(gpsList);
    frameDurs = durList(ndx);

return