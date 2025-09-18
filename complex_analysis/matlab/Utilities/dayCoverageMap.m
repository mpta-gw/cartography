function [n]=dayCoverageMap(segments,siderealtimes)
%  function [n]=dayCoverageMap(segments,times)
%
%  dayCoverageMap  -- returns the numer of segments that cover 
%                     each sidereal time in vector siderealtimes
%                    
%  arguments: segmens: list of segments, format: [GPS_starttime, duration]
%             siderealtimes: vector of sidereal times
%  output:    n: nomber of segments that cover each sidereal time
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=size(segments);
if s(2) == 2
  Nseg=s(1);
else
  error('1st argument must by Nx2 array with a list of GPS start times and durations');
end

% create n with right dimension
n=zeros(size(siderealtimes));

% get sidereal start and stop times, don't apply the mod yet
starts=GPStoGreenwichMeanSiderealTime(segments(:,1),false)/24;
ends=  GPStoGreenwichMeanSiderealTime(segments(:,1)+segments(:,2),false)/24;
durs=ends-starts;

% take care of segments that are >=1sidereal day, increment segment counter for all times
fulldays=floor(durs);
durs=durs-fulldays;
n=n+sum(fulldays);

% take the modulus of start, end and requested times
starts=mod(starts,1);
ends  =mod(ends,1);
times =mod(siderealtimes/24,1);

% calculate for each time whether is inside or outside the segment
% idea: there are 3 relevant times: start (ts), end (te), and requested time (t)
% there are 6 possible permutations:
% 
% time order  ts<te  ts<t  t<te   t in segment = xor(ts<te,ts<t,t<te)
% -------------------------------------------------------------------
%  ts<te<t     1      1      0      0
%  ts<t <te    1      1      1      1
%  t <ts<te    1      0      1      0
%  te<ts<t     0      1      0      1
%  te<t <ts    0      0      0      0
%  t <te<ts    0      0      1      1
%
% special case: make sure that ts!=te (dur=multiple of full day, e.g. 0sec)
%   (if that is the case: disregard segment - full days have already been counted) 
%   provided ts!=te: this algorithm then excludes the edges of segments, as desired.
%
Lse=starts<ends;
LseSpecial=starts~=ends;
for kk=1:length(times)
  Lst=starts    < times(kk);
  Lte=times(kk) < ends;
  dn=sum(and(xor(xor(Lse,Lst),Lte),LseSpecial));
  n(kk)=n(kk)+dn;
end
