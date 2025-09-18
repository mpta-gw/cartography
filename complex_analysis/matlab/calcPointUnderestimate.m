function [map,dpix]=calcPointUnderestimate(Gamma,res)

Lmax=sqrt(size(Gamma,1))-1;

try,
  RA=360/res;
  DECL=180/res+1;
catch
  res=90/Lmax;
  RA=4*Lmax;
  DECL=2*Lmax+1;
end;

dpix=real(diagPixel(Gamma,res))*4*pi/(Lmax+1)^2;

map=map2plm(dpix,RA,DECL,Lmax);
