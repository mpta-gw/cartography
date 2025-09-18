function [sdl,l,bta,lfit,logy,logpm]=plm2spec(lmcosi,norma)
% [sdl,l,bta,lfit,logy,logpm]=PLM2SPEC(lmcosi,norma)
%
% Calculates the power spectrum of real spherical harmonic
% sine and cosine coefficients contained in the matrix 'lmcosi'
%
% INPUT:
%
% lmcosi               Spherical harmonic coefficients [l m Ccos Csin]
% norma                [1] multiplication by (l+1) 
%                          (Lowes, JGR 71(8), 2179 [1966])
%                          This gives the mean-square value of the
%                          gradient of a potential in Schmidt-harmonics
%                      [2] division by (2*l+1) [default]
%                          (Nagata, JGeomagGeoel 17, 153-155 [1965])
%                      [3] none
%
%
% OUTPUT:
%
% sdl                  Spectral density: energy per degree
% l                    Degree
% bta                  Spectral slope of loglog(l,sdl)
% lfit,logy            Spectral line plot given by loglog(lfit,logy)
% logpm                Error on spectral line plot given by
%                      loglog(lfit,logpm)
%
% EXAMPLE:
%
% [sdl,l,bta,lfit,logy,logpm]=plm2spec(fralmanac('EGM96'));
%
% See also ACTSPEC
%
% The normalization by (2l+1) is what's required when the spherical
% harmonics are normalized to 4pi. See DT p. 858. A "delta"-function then
% retains a flat spectrum. 
% See papers by Hipkin 2001, Kaula 1967, Lowes 1974. 
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 19th, 2004

defval('norma',2)
lmin=lmcosi(1);
lmax=lmcosi(end,1);

pin=0;
for l=lmin:lmax
     clm=shcos(lmcosi,l);
     slm=shsin(lmcosi,l);
     pin=pin+1;
     sdl(pin)=clm(:)'*clm(:)+slm(:)'*slm(:);
end

switch norma
 case 1 
  normfac=(lmin:lmax)+1;
 case 2
  normfac=1./(2*(lmin:lmax)+1);
 case 3 
    normfac=1;
    disp('Not further normalized')
 otherwise
  error('No valid normalization specified')
end

% disp(sprintf('Normalization %i',norma))

sdl=normfac.*sdl;
sdl=sdl(:);
l=lmin:lmax; 
l=l(:);
if lmin==0
  in=3;
elseif lmin==1
  in=2;
else
  in=1;
end
lfit=l(in:end);

if nargout>=3
  % Calculate spectral slope
  [bt,E]=polyfit(log10(lfit),log10(sdl(in:end)),1);
  bta=bt(1);
  [logy,loge]=polyval(bt,log10(lfit),E);
  logy=10.^logy;
  logpm=[logy./(10.^loge) logy.*(10.^loge)];
end
