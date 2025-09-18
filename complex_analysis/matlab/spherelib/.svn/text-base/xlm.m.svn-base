function [X,theta,norms]=xlm(l,m,theta,check,tol)
% [X,theta,norms]=XLM(l,m,theta,check,tol)
%
% Calculates normalized (associate) Legendre functions, DT (B.58/B.60).
%
% INPUT:
%
% l      degree (0 <= l <= infinity) [default: random]
% m      order (-l <= m <= l)        [default: all orders 0<=l]
%        l and m can be vectors, but not both at the same time
% theta  argument (0 <= theta <= pi)    [default: 181 linearly spaced]
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check
% tol    Tolerance for optional normalization checking
%
% OUTPUT:
%
% X      The associated normalized Legendre function at the desired argument(s):
%           as a scalar or a row vector with length(theta) columns, OR
%           as a matrix with length(m) rows and length(theta) columns, OR 
%           as a matrix with length(l) rows and length(theta) columns
% mu     The argument(s), which you might or not have specified
% norms  The normalization matrix, which should be the identity matrix
%
% EXAMPLES:
%
% plot(xlm([0:5],0)')
% plot(xlm(5,[])')
%
% Last modified by fjsimons-at-alum.mit.edu, 13.09.2005

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('theta',linspace(0,pi,181))
defval('check',0)
defval('tol',1e-10)

% Revert back to cos(theta)
mu=cos(theta);

% Error handling common to PLM, XLM, YLM
[l,m,mu,check,tol]=pxyerh(l,m,mu,check,tol);

switch check
  case 0
   % Calculation for m>0
   if prod(size(l))==1 & prod(size(m))==1
     X=(-1)^(-m)*sqrt(2*l+1)/sqrt(2-(m==0))/sqrt(4*pi)*...
       rindeks(legendre(l,mu,'sch'),abs(m)+1);
     if m < 0 
       X=(-1)^(-m)*X;
     end
   elseif prod(size(l))==1 
     X=sqrt(2*l+1)/sqrt(4*pi)*...
       (rindeks(legendre(l,mu,'sch'),abs(m)+1)'*...
	diag(((-1).^(-m))./sqrt(2-(m==0))))';
     for index=find(m<0)
       X(index,:)=(-1)^(-m(index))*X(index,:);
     end
   elseif prod(size(m))==1 
     for index=1:length(l)
       X(index,:)=(-1)^(-m)*sqrt(2*l(index)+1)/sqrt(2-(m==0))/sqrt(4*pi)*...
	   rindeks(legendre(l(index),mu,'sch'),abs(m)+1);
       if m<0
	 X(index,:)=(-1)^(-m)*X(index,:);
       end
     end
   else
     error('Specify valid option')
   end
   norms=[];
 case 1
  % Check normalization 
  norms=pxynrm(l,m,tol,'X');
  % Still also need to get you the answer at the desired argument
  X=xlm(l,m,theta,0);
end


