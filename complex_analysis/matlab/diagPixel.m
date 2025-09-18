function [d,RA,DECL,dOmg,U]=diagPixel(M,res)

% same as covarSpH2Pixel, but only returns the diagonal pixels
% M is the Fisher matrix; d is an array contain the diagonal elements of
% Fisher in the pixel basis.
% U converts X_SpH to X_pixel.

global PIXELCONVERSION;

Lmax=sqrt(size(M,1))-1;
try, res; catch res=90/Lmax; end;
checkPixelConversion(Lmax,res);

NN=size(PIXELCONVERSION.U,1);
d=zeros(NN,1);
for jj=1:NN,
  d(jj)=PIXELCONVERSION.U(jj,:)*M*PIXELCONVERSION.U(jj,:)';
end

if nargout>1
    RA  =PIXELCONVERSION.RA;
    DECL=PIXELCONVERSION.DECL;
    dOmg=PIXELCONVERSION.dOmg;
    U   =PIXELCONVERSION.U;
end
