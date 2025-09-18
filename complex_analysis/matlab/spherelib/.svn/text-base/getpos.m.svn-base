function pos=getpos(varargin)
% pos=GETPOS(ah)
% pos=GETPOS(ah,ind)
%
%
% Gets axis handle position or just
% one of its four elements.
% If input is 1xN vector returns Nx4 or Nx3 positions.
%
% See also GETTIT

defval('ah',gca)

pos=get(varargin{1},'Position');

if iscell(pos)
  vr=length(pos{1});
  pos=reshape([pos{:}],vr',length([pos{:}])/vr')';
end

if nargin>1
  pos=pos(:,varargin{2});
end

