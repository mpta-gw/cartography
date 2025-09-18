function x0=projectPositive(x0)

N=length(x0)/2;

Lmax=sqrt(N)-1;

res=90/Lmax;
for ii=1:3
[map,ra,decl,dOmg]=makemap(x0(1:N)+i*x0(N+1:2*N),res);
%max(map)
plm=map2plm(max(map,0),ra,decl,Lmax);
%norm(plm)
x0=[real(plm);imag(plm)];
end
