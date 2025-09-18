function imagef(c11,cmn,matrix)
% IMAGEF(c11,cmn,matrix)
%
% Is the way image and imagesc should've been written.
% FOR GEOGRAPHICAL COORDINATES!
%
% Give matrix file, and the coordinates of element
% matrix(1,1) and matrix(m,n); 
% PIXEL CENTERED COORDINATES
% Therefore, should C11 and CMN give the boundaries
% of the grid, use IMAGEF(c11+[dlo -dla]/2,cmn+[-dlo dla]/2,both)
% with DLO and DLA the X and Y interval
%
% See IMAGEFDIR
%
% Last modified by fjsimons-at-alum.mit.edu, September 17th, 2003

defval('c11',[0 90])
defval('cmn',[360 -90])

cm1=[c11(1) cmn(2)];
c1n=[cmn(1) c11(2)];
imagesc([cm1(1) c1n(1)],[cm1(2) c1n(2)],flipud(matrix));  axis xy

% Equivalent pcolor
% pcolor(linspace(c11(1),cmn(1),size(both,2)+1),...
%    linspace(c11(2),cmn(2),size(both,1)+1),adrc(matrix));

if all(c11==[0 90]) & all(cmn==[360 -90])
  set(gca,'ytick',[-90:30:90])
  set(gca,'xtick',[0:90:360])
  deggies(gca)
  disp('IMAGEF: We made the decision for you')
end

