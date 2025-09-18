function figdisp(name,ext,opt)
% FIGDISP(name,ext,opt)
%
% Suggests print command for figures.
% Assumes environment variable 'EPS' is set.
% Use PAINTERS to override, sometimes.
%
% Last modified by fjsimons-at-alum.mit.edu, June 10th 2004

[p,n]=star69;

defval('name',n)
defval('ext',[])
defval('opt',[])

if ~isstr(ext)
  ext=num2str(ext);
end

if ~isempty(ext)
  fname=fullfile(getenv('EPS'),sprintf('%s_%s',name,ext));
else
  fname=fullfile(getenv('EPS'),sprintf('%s',name));
end

if ~isempty(opt)
  disp(sprintf('print(''-depsc'',''%s'',''%s'')',opt,fname))
else
  disp(sprintf('print(''-depsc'',''%s'')',fname))
end

