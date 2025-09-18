function d=WignerDMatrix(l,beta);

% calculates the rotation matrix for spherical harmonics
% 
% beta is in radians

m=-l:l;
mp=m';


N=(rfactorial(l+mp).*rfactorial(l-mp))*(rfactorial(l+m).*rfactorial(l-m));

S=zeros(size(N));

for s=0:2*l
  mpms = mp*ones(1,2*l+1) - ones(2*l+1,1)*m + s;
  prefac=((invfact(s)*invfact(l-mp-s))*invfact(l+m-s)).*invfact(mpms).*(-1).^mpms;
  CC = cos(beta/2).^(2*l-mpms-s);
  SS = sin(beta/2).^abs((mpms+s));
  S=S+prefac.*CC.*SS;
end

d=sqrt(N).*S;


return

function invfactn = invfact(n)
  invfactn=1./rfactorial(n+(n<0)*1e100);
return

% version of matlab factorial function that keeps y the same shape as x
function y = rfactorial(x)
  y = reshape(factorial(x),size(x));
return
