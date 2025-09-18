function [phix,phi,wx,w,philim,wlim]=graphs(tipe,nvm)
% [phix,phi,wx,w,philim,wlim]=GRAPHS(tipe,nvm)
%
% Computes graphs of wavelets and scaling functions.
% To get dual, switch (f0,f1) and (h0,h1)
%
% See Mallat (1998), p258
%
% Last modified by fjsimons-at-alum.mit.edu, March 31st, 2003

defval('tipe','CDF')
defval('nvm',[2 4])

[h0,f0]=wc(tipe,nvm);

% Number of vanishing moments is equal to (l+1)/2
[h0,h1,f0,f1,l]=prodco(f0,h0);

% Number of iterations
j=7;
% Number of points warranted
N=2^(j+5);
% Calculate number of coefficients
[an,dn,a,d]=dnums(N,h0,f0,j);
% Specify shift (not on the edge)
% Needs to be somewhere in the middle
% Influences the position of the return
n=ceil(an(end)/2);

% Centered unstretched unshifted time axis
tax=linspace(-floor((N+1)/2)+1,N-floor((N+1)/2),N)/2^j;

% Compute scaling function by inverse wavelet transform
a(n)=1;
x=iwt(a,d,an,dn,tipe,nvm);
x=[x{:}];
phi=sum(x,2);
% Figure out the finite support
supp=find(abs(phi)>1e-10);

% Time axis where supported
phix=indeks(tax,supp);
% Scaling function corrected to zeroth level ('father')
phi=indeks(phi,supp)*sqrt(2^j);
% Theoretical support
philim=[0 l];

[an,dn,a,d]=dnums(N,h0,f0,j);

% Compute wavelet at scale 2^j and shift n
d{j}(n)=1;
x=iwt(zeros(1,an(end)),d,an,dn,tipe,nvm);
x=[x{:}];
w=sum(x,2);
% Figure out the finite support
supp=find(abs(w)>1e-10);

% Time axis where supported
wx=indeks(tax,supp);
% Wavelet corrected to zeroth level ('mother')
w=indeks(w,supp)*sqrt(2^j);
% Theoretical support
wlim=[(1-l)/2 (l+1)/2];

