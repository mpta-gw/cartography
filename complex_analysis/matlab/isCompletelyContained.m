function out=isCompletelyContained(GPS,duration,list)
%function out=isCompletelyContained(GPS,list)
%
% Checks whether GPS times (+duration(s)) is totally contained in list
%
% input:    GPS:      Array of GPS times
%           duration: Either one duration or Array of durations
%           list      Segment list of format SegNum, Start, Stop Dur
%                     only Start and Stop are used
%
% output:   out:      boolean array the size of GPS
%
%  Routine written by Stefan Ballmer
%  Contact sballmer@ligo.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


out=zeros(size(GPS));
N=size(list,1);

GPSe=GPS+duration;
for ii=1:N
  out=or(out,and(list(ii,2)<=GPS,GPSe<=list(ii,3)));
end
