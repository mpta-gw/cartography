function [covarP,RA,DECL,dOmg,U]= covarSpH2Pixel(covar,deg)

% transform a covariance matrix from spherical harmonisc basis
% to pixel basis



global PIXELCONVERSION;

Lmax=sqrt(size(covar,1))-1;
try, deg; catch deg=90/Lmax; end;
checkPixelConversion(Lmax,deg);


covarP=PIXELCONVERSION.U*covar*PIXELCONVERSION.U';

if nargout>1
    RA  =PIXELCONVERSION.RA;
    DECL=PIXELCONVERSION.DECL;
    dOmg=PIXELCONVERSION.dOmg;
    U   =PIXELCONVERSION.U;
end
