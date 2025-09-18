function [sigma,sigmaPix,RA,DECL,dOmg,U]=getSigmaMap(covarM,res,RA,DECL,dOmg,U)
% return the map corresponing to the standard deviation for each pixel
% This map is defined as the sqrt(diag(covariancematrix)),
% in the pixel basis
% 
%

Lmax=sqrt(size(covarM,1))-1;
try, res; catch res=90/Lmax; end;

if exist('U')
  NN=size(U,1);
  diagCovarMPix=zeros(NN,1);
  for jj=1:NN,
    diagCovarMPix(jj)=U(jj,:)*covarM*U(jj,:)';
  end
else
    if nargout>2
        [diagCovarMPix,RA,DECL,dOmg,U]=diagPixel(covarM,res);
    else
         diagCovarMPix                =diagPixel(covarM,res);
    end
end
sigmaPix=sqrt(real(diagCovarMPix));
sigma=map2plm(sigmaPix,ceil(360/res),ceil(180/res)+1,Lmax);
