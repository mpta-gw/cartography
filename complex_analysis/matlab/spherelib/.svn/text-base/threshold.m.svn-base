function [c,cnz]=threshold(d,dn,kind,div)
% [c,cnz]=threshold(d,dn,kind,div)
%
% Thresholds wavelet or scaling coefficients
% by selection (see Mallat p. 437 ff)
%
% INPUT:
%
% d       -  Coefficient cell array
% dn      -  Length of the cell array
% kind    -  'soft' or 'hard'
%            (May need to look into halving T for soft thresholding.)
% div     -  Divides threshold level, if at all needed
%            Higher values imply lower threshold
%
% OUTPUT:
%
% c       -  Thresholded coefficients
% cnz     -  Number of nonzero coefficients in every cell
%
% Last modified by fjsimons-at-alum.mit.edu, December 04th 2003

defval('kind','soft')
defval('div',1)

for index=1:length(dn)
  % Calculate location estimate
  % Calculate scale estimate for coefficients at all resolution scales
  % Use median absolute deviation from the median
  % Mallat Equation 10.53
  dloc(index)=median(d{index});
  dscal(index)=median(abs(d{index}-dloc(index)))/0.6745;
  T(index)=dscal(index)*sqrt(2*log(dn(index)))/div;
  switch kind
   case 'hard'
    dist=abs(d{index})>T(index);
    c{index}=d{index}.*dist;
    if index==1 ; disp('Hard thresholding') ; end
   case 'soft'
    dist=(abs(d{index})-T(index));
    dist=(dist+abs(dist))/2;
    c{index}=sign(d{index}).*dist;
    if index==1 ; disp('Soft thresholding') ; end
  end
  % Return effective number of (nonzero) coefficients
  cnz(index)=sum(~~dist);
end





