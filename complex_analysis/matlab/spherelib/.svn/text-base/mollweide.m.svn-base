function [x,y]=mollweide(lon,lat,lonref)
% [x,y]=MOLLWEIDE(lon,lat,lonref)
%
% X and Y are the coordinates of LON and LAT
% after Mollweide projection with respect to
% LONREF. Everything in radians.
%
% Last modified by fjsimons-at-alum.mit.edu , June 27th, 2003

defval('lonref',pi)

i=(lonref-lon)>pi;

lon(i)=lon(i)+2*pi;
i=(lon-lonref)>pi;
lon(i)=lon(i)-2*pi;

theta=lat;

dtheta=1;
while max(abs(dtheta))>1e-8,
	dtheta=-(theta+sin(theta)-pi*sin(lat))./(1+cos(theta));
	theta=theta+dtheta;
end
theta=theta/2;

x=2*sqrt(2)*(lon-lonref).*cos(theta)/pi;
y=sqrt(2)*sin(theta);
