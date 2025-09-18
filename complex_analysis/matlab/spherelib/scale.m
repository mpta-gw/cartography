function scalmat=scale(inmat,ranscale)
% scalmat=SCALE(inmat,ranscale)
% 
% Rescales a matrix so its range is now 'ranscale'
% As an example,:
% plot(scale(scale(sin(0:0.01:3*pi),[1898 9100]),[-1 1]),'r-')
%
% Last modified by fjsimons-at-alum.mit.edu, August 12th 1998

if diff(ranscale)<=0
  ranscale=fliplr(ranscale);
end

reensj=max(inmat(:))-min(inmat(:));
minim=min(inmat(:));

scalmat=ranscale(1)+(inmat-minim)*range(ranscale)/reensj;

