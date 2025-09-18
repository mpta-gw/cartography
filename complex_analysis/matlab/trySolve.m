% This is a demo script that was used to explore possible
% ways of regularizing the deconvolution step
% both with and without assuming positivity


global NUMTOOLS
global RAWDATA

p=load('test_SpH_SkyMap_24jobs_short_ra6_dec45_H2e46_L10_time.mat');
%p=load('test_SpH_SkyMap_24jobs_short_monopole_H2e46_L10_time.mat');
Lmax=10;
p.x=subX(p.x./1e46,Lmax);
p.fisher=subFisher(p.fisher./(1e46)^2,Lmax);
p.invfisher=inv(p.fisher);
p.p=p.invfisher*p.x;

RAWDATA=p;


% regularization
ii=floor(80/121*(Lmax+1)^2);
[U,S,V] =svd(p.fisher);
eg=diag(S);
eg(ii+1:end)=eg(ii);
RAWDATA.fisher=U*diag(eg)*V';
RAWDATA.invfisher=V*diag(1./eg)*U';

Pcovar=V*diag(1./eg)*S*diag(1./eg)*U';


eta=diag(1./eg)*U'*p.x;
sigmaeta=sqrt(diag(S))./eg;


N=length(p.x);
lmax=sqrt(N)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize with best solution, ignoring positivity
data=inv(RAWDATA.fisher)*p.x;
x0=[real(data);imag(data)];
NUMTOOLS.BESTVALUE=minimizeThis(x0);

% project this to positive
x0=projectPositive(x0);

%alternative initialization
%data=zeros(size(p.x));
%x0=[real(data);imag(data)];


NUMTOOLS.func     ='minimizeThis';
NUMTOOLS.stepsize = 1e-3*ones(size(x0));
NUMTOOLS.maxiter  = 100;
NUMTOOLS.maxerr   = 1e-60;

x0=projConjGradSolve(x0);

d=x0(1:N)+i*x0(N+1:2*N);

datainv=inv(p.fisher)*p.x;
datainv1=inv(RAWDATA.fisher)*p.x;

plotPlmMap(p.x);
plotPlmMap(datainv);
plotPlmMap(datainv1);
title('Best fit')
plotPlmMap(d);
title('Best positive fit')

[sigma,sigmaPix]=getSigmaMap(Pcovar,1);
plotMapAitoff(sigmaPix,360,181,-1);
title('Standard Deviation')
[dPix    ,RA,DECL,dOmg] = makemap(d       ,1);
plotMapAitoff(dPix./sigmaPix,360,181,-1);
title('Best positive fit SNR');

