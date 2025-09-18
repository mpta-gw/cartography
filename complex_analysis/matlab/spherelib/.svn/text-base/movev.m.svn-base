function varargout=movev(fax,amt)
% MOVEV(fax,amt)
% [new,old]=MOVEV(...)
%
% Moves an axis handle vertically by an 'amt'.
% Also works for non-axis objects.
% + is up, - is down
%
% INPUT:
%
% fax      Axis handle
% amt      Amount of move
%
% Returns new and old positions.
%
% Last modified by fjsimons-at-alum.mit.edu, 25.2.2005

defval('amt',0.1)

if nargout>=1
  varargout{1}=getpos(fax);
end

for index=1:length(fax)
  oldpos=get(fax(index),'Position');
  newpos=repmat(0,1,length(oldpos));
  newpos(2)=amt;
  newpos=oldpos+newpos;
  set(fax(index),'Position',newpos)
end

if nargout>=2
  varargout{2}=getpos(fax);
end
