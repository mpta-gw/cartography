function [covar]= covarPixel2SpH(covarP,RA,DECL,Lmax)

% transform a covariance matrix from 
% from pixel basis to spherical harmonisc basis

if 1
  plmHalf=map2plm(covarP  ,RA,DECL,Lmax,1);
  plm    =map2plm(plmHalf',RA,DECL,Lmax,1);

  covar=plm';
else
  % alternative code; works, but has less numerical precision due to crude integration
  [lvec,mvec]=getLMvec(Lmax);

  N=length(lvec);

  %try, deg; catch deg=90/Lmax; end;
  deg=360/RA;

  for ii=1:N
    ylm=zeros(N,1); ylm(ii)=1;
    [map,RA,DECL,dOmg] = makemap(ylm,deg,0,1);
    U(1:length(map),ii)=map.*dOmg;
  end
  covar=U'*covarP*U;
end

