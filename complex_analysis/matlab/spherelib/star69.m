function [path,name,ext,ver]=star69
% [path,name,ext,ver]=STAR69
%
% Returns calling function name if exists.
%
% Last modified by fjsimons-at-alum.mit.edu, 10.11.2004


v=version;

str=dbstack;

if v(1)== '5';
  if length(str) > 2
    str=str(3).name;
  else
    str='';
  end
elseif v(1)== '6'
  str=str(length(str)).name;
elseif v(1)== '7'
  str=str(length(str)).name;
end

[path,name,ext,ver]=fileparts(str);

