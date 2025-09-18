function [Nm,Nsum]=nsubm(N,m,method,L)
% [Nm,Nsum]=NSUBM(N,m,method,L)
%
% Calculates Nm, the partial Shannon number at a given
% angular order m, from the full Shannon number, N, by a variety
% of (approximate) formulations. Returns the sum over all angular
% orders. 
%
% Use this to figure out how many tapers to take at a particular angular
% order if you have the full Shannon number (the area of the region and
% the bandwidth will calculate this for you). This is useful for the
% axisymmetric single-order cases where you'll want to calculate the
% tapers using Grunbaum's fast method and you might not bother with the
% exact eigenvalues at all.
%
% INPUT:
%
% N         Full Shannon number
% m         Maximum angular order; all will be computed
% method    1 Using Bessel functions (asymptotic)
%           2 Using hypergeometric functions (asymptotic; slower)
%           3 Exact formula using GL integration (requires L)
% L         Bandwidth (only required for exact formulation)
%
% OUTPUT:
%
% Nm        The partial Shannon number at angular order m
% Nsum      The sum of all Nm taking into account degeneracy
%
% EXAMPLE:
%
% L=18; TH=20; N=(L+1)^2*(1-cos(TH/180*pi))/2;
% [Nm,Nsum]=nsubm(N,L,3,L); N-Nsum
%
% Last modified by fjsimons-at-alum.mit.edu, 23.2.2005

defval('N',30)
defval('m',30)
defval('method',1)
defval('L',m)

switch method
  case 1
   s=2*sqrt(N);
   J=besselj(0:m+1,s);
   Nm=2*N*[J(1:m+1).^2+J(2:m+2).^2]-...
      (2*[0:m]+1)*sqrt(N).*J(1:m+1).*J(2:m+2)-...
      1/2*[0:m].*[1-J(1)^2-2*cumsum([0 J(2:m+1).^2])];
   %disp('Asymptotic (Bessel) formalism')
 case 2
  if(~isempty(help('symbolic')))
      for ml=0:m
	Nm(ml+1)=N^(ml+1)*...
		  double(hypergeom([1/2+ml,1+ml],[1+2*ml,2+ml,2+ml],-4*N))/...
		  (1+ml)/gamma(2+ml)/gamma(1+ml);
      end
      disp('Asymptotic (hypergeometric) formalism')
    else
      error('Need symbolic math toolbox')
    end
 case 3
  % The area, must be smaller than 4 pi of course
  A=4*pi*N/(L+1)^2;
  % What's the rough colatitude
  th0=acos(1-A/2/pi);
  if ~isreal(th0)
    Nm=repmat(NaN,1,m+1);
    warning('Area must be smaller than unit sphere; NaNs returned')
  else
    Nm=repmat(0,1,m+1);
    if L<m
      warning('Degree must be bigger or equal than order; zeroes returned')
    end
    for m1=0:m
      for L1=m1:L
	% Need Schmidt*(2L+1)/4/pi for fully normalized, but want Xlm, not Ylm
	Nm(m1+1)=Nm(m1+1)+legendreprodint(L1,m1,L1,m1,cos(th0),'gl')...
		 *(2*L1+1)/4/pi/(2-(m1==0));
      end
    end
    Nm=Nm*2*pi;
    disp(sprintf('Gauss-Legendre formalism with th0= %3.3i',...
    		 round(th0*180/pi)))
  end
 otherwise
  error('Specify valid option')
end

Nsum=Nm(1)+2*sum(Nm(2:end));
