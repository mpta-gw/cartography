function   mu_UL=getUpperLimit(Y,sigma,confidence,calError)
%function   mu_UL=getUpperLimit(Y,sigma,confidence,calError)
%
% marginalizes over the calibration uncertainty assuming
% a normal distribution and calculates a Bayesian upper limit
%
% Y         : point esimate (vector)
% sigma     : theoretical sigma (vector)
% confidence: confidence level for upper limit (e.g. 0.9)
% calError  : Calibration uncertainty
%              (for 2 interferometers: =sqrt(err1^2+err2^2) )
%
% Author: Stefan Ballmer, sballmer@ligo.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(Y);
mu_UL=zeros(size(Y));
for n=1:N
  Omg=sigma(n)*(0:1e-3:10);
  signew2=(sigma(n).^2+calError^2*Omg.^2);
  p=sqrt(1./signew2).*exp(-1/2*(Y(n)-Omg).^2./signew2);
  % this is equivalent to (up to overall scaling)
  %  for k=1:length(Omg)
  %    F=@(l) exp(-1/2*(Y-(1+l)*Omg(k)).^2/sigma^2).*exp(-1/2*l.^2/calError^2);
  %    p(k)=quad(F,-100*calError,100*calError);
  %  end
  cp=cumsum(p,2);
  cp=cp./cp(end);
  in2=find(cp<0.999);
  mu_UL(n)=interp1(cp(in2),Omg(in2),confidence);
  if mod(n,500)==0
    fprintf('n= %d\t UL(f)= %g\n',n,mu_UL(n));
  end
end
%toc
