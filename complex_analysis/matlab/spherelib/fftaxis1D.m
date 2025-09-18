function [fax,selekt]=fftaxis1D(sig,nfft,spun)
% [fax,selekt]=fftaxis1D(sig,nfft,spun)
%
% See also FFTAXIS
%
% For 1-D signals, returns a selection vector 'selekt'
% and a frequency axis 'fax' based on the total dimension-having
% length of the signal 'spun' (as in: 123 m)
% So axis for ISREAL is from 0 to Nyquist, for ISCOMPLEX from 0 to twice
% that (as an alternative to from -f_N to f_N)
%
% Use without FFTSHIFT, as in 
% FX=fft(S,nfft); SX=abs(FX).^2;
% plot(fax,SX(selekt)

% Last modified by fjsimons-at-alum.mit.edu, May 8th, 2004

xsint=spun/(length(sig)-1);

if isreal(sig)
  if rem(nfft,2); selekt = [1:(nfft+1)/2]; else; selekt = [1:nfft/2+1]; end
else  
  selekt = 1:nfft; 
end
fax = (selekt - 1)'/xsint/nfft;
