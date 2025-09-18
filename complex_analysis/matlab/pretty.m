function pretty(varargin)
% function pretty()
% E. Thrane
% makes plots look good/readable for publications
% the font sizes are set to a standard size
if (nargin == 1)
  fs = varargin{1};
else
  fs = 20;
end

h_xlabel = get(gca, 'XLabel');
set(h_xlabel, 'FontSize', fs);
h_ylabel = get(gca, 'YLabel');
set(h_ylabel, 'FontSize', fs);
h_zlabel = get(gca, 'ZLabel');
set(h_zlabel, 'FontSize', fs);
h_title = get(gca, 'Title');
set(h_title, 'FontSize', fs);
set(gca, 'FontSize', fs);

return;
