function plm=randPlm(Lmax,doNorm)

% function returns a random realisation of spherical harmonics
% if doNorm is true (default), the vector will be normalised
% to a norm of 1, otherwise the vector has a norm of 1 on average.

try, doNorm; catch doNorm=true; end;

N=(Lmax+1)*(Lmax+2)/2;

Nindep=(Lmax+1)^2;


plm=plmreal2plm(randn(N,2)./(Lmax+1));

if doNorm
  plm=plm./norm(plm);
end
