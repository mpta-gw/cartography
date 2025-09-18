function lmrlim=rsh2cpx(lmcosi)
% lmrlim=RSH2CPX(lmcosi)
%
% Transforms REAL spherical harmonics into COMPLEX ones
%
% INPUT:
%
% lmcosi        Standard matrix listing l, m, Ccos and Csin
%
% OUTPUT:
%
% lmrlim        Matrix listing coefficients for m=-l:l as
%               l abs(m) Creal Cimag
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 12, 2004

[l,m,mz]=addmon(max(lmcosi(:,1)));

lmrlim=lmcosi;

% Divide everything by sqrt(2)...
lmrlim(:,3:4)=lmcosi(:,3:4)/sqrt(2);
% ... except the m=0 component
lmrlim(mz,3:4)=lmcosi(mz,3:4);




