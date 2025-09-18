function [xf,co,npol,npas,tipe]=lowpass(x,Fs,co,npol,npas,tipe)
% [xf,co,npol,npas,tipe]=LOWPASS(x,Fs,co,npol,npas,tipe)
%
% Filters signal 'x' with filter 'tipe' and corner
% frequency 'co' in Hz with 'npol' number of poles and in 
% 'npas' passes. Sampling rate is 'Fs' in Hz.
%
% Compare in SAC lp butter co 5 n 2 p 1
%
% INPUT:
%
% x         The signal
% Fs        Its sampling frequency
% co        The corner frequency
% npol      The number of poles
% npas      The number of passes
% tipe      The filter name
% trending  'linear' or 'constant'
%
% OUTPUT:
%
% xf      The filtered signal
%
% Last modified by fjsimons-at-alum.mit.edu, 22.11.2004

defval('Fs',110)
defval('co',5)
defval('npol',2)
defval('npas',1)
defval('tipe','butter')
defval('trending','linear')

% Corner frequency is in Hertz, now it is as a fraction of
% half the sampling rate.
Wn=2*co/Fs;
[B,A]=feval(tipe,npol,Wn);

xf=filter(B,A,detrend(x(:),trending));

if npas==2
  xf=flipud(filter(B,A,detrend(flipud(xf(:)),trending)));  
end

disp(sprintf('LOWPASS %3.3f Hz %i pass %i poles %s',co,npas,npol,tipe))
