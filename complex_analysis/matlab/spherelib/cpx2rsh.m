function lmcosi=cpx2rsh(lmrlim)
% lmcosi=CPX2RSH(lmrlim)
%
% Transforms COMPLEX spherical harmonics into REAL ones
%
% INPUT:
%
% lmrlim        Matrix listing coefficients for m=-l:l as
%               l abs(m) Creal Cimag
%
% OUTPUT:
%
% lmcosi        Standard matrix listing l, m, Ccos and Csin
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 12, 2004

[l,m,mz]=addmon(max(lmrlim(:,1)));

lmcosi=lmrlim;

% Multiply everything by sqrt(2)...
lmcosi(:,3:4)=lmrlim(:,3:4)*sqrt(2);
% ... except the m=0 component
lmcosi(mz,3:4)=lmrlim(mz,3:4);




