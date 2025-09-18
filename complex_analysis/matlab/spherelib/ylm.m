function [Y,theta,phi,norms]=ylm(l,m,theta,phi,check,tol)
% [Y,theta,phi,norms]=YLM(l,m,theta,phi,check,tolb)
%
% Calculates normalized real spherical harmonics, DT (B.72).
%
% INPUT:
%
% l      degree (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders 0<=l]
%        l and m can be vectors, but not both at the same time
% theta  colatitude (0 <= theta <= pi)   [default: 181 linearly spaced]
% phi    longitude (0 <= theta <= pi)    [default: 361 linearly spaced]
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check
% tol    Tolerance for optional normalization checking
%
% OUTPUT:
%
% Y      The real spherical harmonics at the desired argument(s):
%           As matrix with dimensions of 
%           length(theta) x length(phi) x max(length(m),length(l))
% theta  The latitude(s), which you might or not have specified
% phi    The longitude(s), which you might or not have specified
% norms  The normalization matrix, which should be the identity matrix
%
% EXAMPLES:
%
% plotplm(ylm(2,-1,[],[],1),[],[],1); colormap(flipud(gray(7)))
%
% Last modified by fjsimons-at-alum.mit.edu, 13.09.2005

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('theta',linspace(0,pi,181))
defval('phi',linspace(0,2*pi,360))
defval('check',0)
defval('tol',1e-10)

% Make sure phi is a row vector
phi=phi(:)';

% Revert back to cos(theta)
mu=cos(theta);

% Error handling common to PLM, XLM, YLM
[l,m,mu,check,tol]=pxyerh(l,m,mu,check,tol);

% Straight to the calculation, check normalization on the XLM
[X,theta,norms]=xlm(l,abs(m),theta,check);

% Initialize the matrix with the spherical harmonics
Y=repmat(NaN,[length(theta) length(phi) max(length(m),length(l))]);

% Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not
P=diag(sqrt(2-(m(:)==0)))*...
  cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi)));

% Make the real spherical harmonics
if prod(size(l))==1 & prod(size(m))==1
 Y=X'*P;
else
  for index=1:max(length(m),length(m))
    Y(:,:,index)=X(index,:)'*P(index,:);
  end
end

