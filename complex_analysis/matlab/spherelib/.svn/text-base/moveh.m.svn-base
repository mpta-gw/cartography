function varargout=moveh(fax,amt)
% MOVEH(fax,amt)
% [new,old]=MOVEH(...)
%
% Moves an axis handle horizontally by an 'amt'.
% Also works for objects.
% + is to the right
% - is to the left
%
% Returns new and old positions.
%
% Modify later to include
% set(fax(index),'Vertices',get(fax(index),...
%     'Vertices')-repmat([1 0 0],4,1))
%
% Last modified by fjsimons-at-alum.mit.edu, 25.2.2005

defval('amt',0.1)

if nargout>=1
  varargout{1}=getpos(fax);
end

for index=1:length(fax)
  oldpos=get(fax(index),'Position');
  newpos=repmat(0,1,length(oldpos));
  newpos(1)=amt;
  newpos=oldpos+newpos;
  set(fax(index),'Position',newpos)
end

if nargout>=2
  varargout{2}=getpos(fax);
end

