function prefix=pref(string,varargin)
% prefix=PREF(string)
% prefix=PREF(string,delimiter)
%
% Finds the prefix or head in a string, say a filename.
% The delimiter is a '.' (period) by default, but may be set.
%
% See also SUF
%
% Last modified by fjsimons-at-alum.mit.edu, 17.11.2004

if nargin==1
  delim= '.';
else
  delim=varargin{1};
end  

prefix=strtok(string,delim);





