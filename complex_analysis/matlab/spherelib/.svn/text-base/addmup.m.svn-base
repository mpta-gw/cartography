function nrmorl=addmup(nr,allertour)
% nrmorl=ADDMUP(nr,direction)
%
% Calculates the number of real spherical harmonic orders
% that belong to an expansion from degree l=0 to L - or vice versa. 
% Since m=0:l for real spherical harmonics, the number of orders per
% degree is l+1. This program thus calculates
% $\sum\limits_{l=0}^{L}(l+1)$ and its inverse, for a given L, using
% analytical expressions. Works for matrices, vectors as well as
% scalars. Note: degrees 0 and 1 are of course included.
%
% INPUT:
%
% nr, 'a' Degree of the expansion; calculates number of orders (default)
% nr, 'r' Number of orders; calculates the degree of the expansion 
%
% OUTPUT:
%
% nrmorl   Depending on 'direction', the number of m's or the degree L
%
% EXAMPLE: 
%
% addmup(3)==[1+2+3+4]
%
% Last modified by fjsimons-at-alum.mit.edu, 22.2.2005

defval('allertour','a')

switch  allertour
  case 'a'
    nrmorl=1/2*(nr+1).^2+1/2*nr+1/2;
  case 'r'
    nrmorl=-3/2+1/2*sqrt(1+8*nr);
    if prod((round(nrmorl)==nrmorl)+0)==0
      warning('Invalid entry - noninteger result')
    end
 otherwise
  error('Specify a valid option')
end


