function newpos=shrink(handle,widthfactor,heightfactor)
% newpos=SHRINK(handle,widthfactor,heightfactor)
%
% This function reduces the width of an object 'handle' by 
% a factor 'widthfactor' and its height by 'heightfactor'
%
% If handle is a vector and width and height scalars,
% applies the same scaling to all the handles.

% Last modified by fjsimons-at-mit.edu, April 6th, 2000

property= 'Position';
handle=handle(:);
if prod(size(widthfactor))*prod(size(widthfactor))~=1
  error('Cannot deal with this multiple input yet')
end
for index=1:length(handle)
  oldvalue=get(handle(index),property);
  newvalue=[oldvalue(1)+oldvalue(3)*((1-1/widthfactor)/2)  ...
	oldvalue(2)+(oldvalue(4)*(1-1/heightfactor))/2 ... % Division by 2 is new
	oldvalue(3)/widthfactor ...
	oldvalue(4)/heightfactor];
  set(handle(index),property,newvalue);
  if nargout>0
    newpos(index,:)=newvalue;
end
end