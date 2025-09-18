function U = pixel2sph(lmax, deg)
%function U = pixel2sph(lmax, deg)
% to convert the covar matrix to sph from pixel:
%   pCovar_sph=U'*pCovar_pixel*U;
%   X_sph = U'*X_pixel;
% by E Thrane based on code by S Ballmer

lmax=strassign(lmax);
deg=strassign(deg);

[lvec,mvec]=getLMvec(lmax);
N=length(lvec);
for ii=1:N
  ylm=zeros(N,1); ylm(ii)=1;
  [map,RA,DECL,dOmg] = makemap(ylm,deg,0,1);
  U(1:length(map),ii)=map.*dOmg;
end

try
  save(['pixel2sph_' num2str(lmax) '_' num2str(deg) '.mat'], 'U');
catch
  save('pixel2sph_abort.mat');
end
