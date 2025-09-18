function [a,d,an,dn,ts,cf]=wt(x,tipe,nvm,n,pph,intel)
% [a,d,an,dn,ts,cf]=wt(x,tipe,nvm,n,pph,intel)
% 
% For input signal 'x' and filters 'h0' and 'f0' 
% (part of a suitable product filter) performs
% a wavelet decomposition by 'n' iterations of filters.
%
% Returns cell arrays with the coefficients.
%
% OUTPUT
%
% 'a'   Approximation/scaling coefficients (after 'n' lowpasses)
% 'd'   Details/wavelet coefficients (after successive highpasses)
% 'an'  Number of approximation coefficients at each level
% 'dn'  Number of detailed coefficients at each level
% 'ts'  How long it took to calculate this (s)
% 'cf'  Compression factor (by the cumulative number of coefficients)
% 
% INPUT
%
% 'x'      The signal to be analyzed
% 'tipe'   'Daubechies' or 'CDF'
% 'nvm'    Number of vanishing (primal & dual) moments 
% 'n'      Number of filter bank iterations
% 'pph'    Method of calculation
%          1= Time-domain full bitrate (inefficient);
%          2= Time-domain polyphase (inefficient);  
%          3= Z-domain polyphase (default)
%          4= Lifting
%
% Last modified by fjsimons-at-alum.mit.edu, 3/18/2002

t0=clock;

% Default is the 5-scale CDF 2/4 construction on Doppler noise
defval('x',real(dopnoise(500,200,60,10,70,128)))
defval('tipe','CDF')
defval('nvm',[2 4])
defval('n',5)
defval('pph',3)
defval('intel',0)

if mod(length(x),2)
  warning('Reconstruction will not work for odd length array')
  pause(10)
end

[h0,f0,P,U,Kp,Ku]=wc(tipe,nvm);

% Need to figure out how to do this recursively
a=x(:);
for index=1:n
  [a,d{index}]=abank(a,h0,f0,P,U,Kp,Ku,pph,intel);
  an(index)=length(a);
end

% Know in advance what the lengths are; except the top level,
% which you'll have to assume is even and given by 2*dn{1}.
dn=cellfun('size',d,1);

ts=etime(clock,t0);
disp(sprintf(' Analysis took %8.4f s',ts))

nuco=cumsum([an(end) fliplr(dn)]);
cf=100*nuco/length(x);
