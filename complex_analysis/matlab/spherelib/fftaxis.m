function varargout=fftaxis(fsize,ftsize,spun)
% [xfaks,yfaks,fnx,fny,xsint,ysint]=FFTAXIS(fsize,ftsize,spun)
% [xfaks,fnx,xsint]=FFTAXIS(fsize,ftsize,spun)
%
% Makes linear frequency axis going through zero.
% -f_N < f <= f_N with f_N=1/(2Dt)
% and floor((dim+1)/2) is at zero.
%
% For 1-D and 2-D but see also FFTAXIS1D
%
% Now works for vectors and matrices with the SIZE command
% fftaxis(size(fx),size(fftx),max(x))
% fftaxis(size(FX),size(FFTX),[lY lX])
%
% 'fsize' is the SIZE of the field (Y X)
% 'ftsize' is the SIZE of its transform (kY kX)
% 'spun' is the physical dimension of the field (lY lX)
%
% A better frequency axis that always goes through zero
% This should return a vector of physical frequencies such that
% -f_N < f <= f_N with f_N=1/(2Dt) the Nyquist rate
% which it returns as 'fnx' and 'fny'
%
% The advantage is that the frequencies are symmetric around the
% zero-center, as they should. The magnitude of the Fourier
% transform of a discrete-time signal is always an even function.
% Analogously, for 2D, the center is a symmetry point.
%
% floor((dim+1)/2) should always give the zero-frequency point
%
% So if we put in
% [a,b,c,d,e,f]=fftaxis([101 51],[100 111],[100 50])
% we see how the frequency goes from -1/2 to 1/2;
% alternatively the angular frequency goes from -pi to pi;
% we can also plot the frequency normalized by the Nyquist:
% in that case, this scale should go from -1 to 1.
%
% The lowest frequency resolvable that is not the DC-component
% is the Rayleigh frequency, given by 1/T=1/NDt with T the data length.
%
% Last verified by fjsimons-at-alum.mit.edu, May 8th, 2004

% From Percival and Walden, 1993

% Figure sampling interval in both directions
if nargout>3
  xsint=spun(2)/(fsize(2)-1);
  ysint=spun(1)/(fsize(1)-1);
  % Calculate centered frequency axis (PW p 112)
  M=ftsize(1);
  N=ftsize(2);
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  intvectY=linspace(-floor((M+1)/2)+1,M-floor((M+1)/2),M);
  xfaks=intvectX/N/xsint;
  yfaks=intvectY/M/ysint;
  % Calculate Nyquist frequencies
  fnx=1/2/xsint;
  fny=1/2/ysint;
  posout={ 'xfaks','yfaks','fnx','fny','xsint','ysint'};
else
  xsint=spun/(max(fsize)-1);

  % Calculate centered frequency axis (PW p 112)
  N=max(ftsize);
  intvectX=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N);
  xfaks=intvectX/N/xsint;
  
  % Calculate Nyquist frequencies
  fnx=1/2/xsint;
  posout={ 'xfaks','fnx','xsint'};
end

for index=1:nargout
  eval([ 'varargout{index}=',posout{index},';'])
end
