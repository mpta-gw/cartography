function [K,kx,ky]=knum2(mn,pl)
% [K,kx,ky]=knum2([M N],[lY lX])
%
% M,N       Size of FFT matrix, equal to size of original matrix
%           [M   N] [rows columns]
% lY,lX     Physical length of original matrix
%           [lY lX] [down across]
%
% Returns wavenumbers belonging to a 2D Fourier
% transform of size (MxN) after FFTSHIFT.
%
% This is only subtly different from FFTAXIS
% The Nyquist rate is (1,1) and the DC rate is
% the "center point" which is in the LR quadrant;
% the way FFTSHIFT expects it.
%
% We used to do this as 
% [xas,yas,fnx,fny,xsint,ysint]=fftaxis(size(ddepth),size(ddepth),[leny lenx]);
% [XAS,YAS]=meshgrid(2*pi*xas,2*pi*yas');
% K=sqrt(XAS.^2+YAS.^2);
% with the problem being that real data sets become complex after being
% manipulated in the wave number domain and inverse transformed. This is
% not the case anymore.
%
% Last modified by fjsimons-at-alum.mit.edu, October 22nd, 2003

M=mn(1);
N=mn(2);

ly=pl(1);
lx=pl(2);

kx=2*pi*linspace(-floor(N/2),floor((N-1)/2),N)/N/lx*(N-1);
ky=2*pi*linspace(-floor(M/2),floor((M-1)/2),M)/M/ly*(M-1);

[KX,KY]=meshgrid(kx(:),ky(:));
K=sqrt(KX.^2+KY.^2);






