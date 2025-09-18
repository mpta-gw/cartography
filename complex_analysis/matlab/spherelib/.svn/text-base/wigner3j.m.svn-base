function W=wigner3j(l1,l2,l3,m1,m2,m3,method)
% W=wigner3j(l1,l2,l3,m1,m2,m3,method)
%
% Calculates absolute values of Wigner 3j coefficients in the context of
% unnormalized associated Legendre functions. One half the inner product
% over the entire interval of three standard, unnormalized, associated
% Legendre polynomials of degree l1,l2,l3 and positive order m1,m2,m3
% equals the square of the corresponding Wigner 3j symbol, 
% as per Dahlen & Tromp Eq. (C.226). 
% 
% See also: http://www-stone.ch.cam.ac.uk/cgi-bin/wigner.cgi
%
% INPUT:
%
% l1,l2,l3     Top row of the Wigner3j symbol; positive angular degree
% m1,m2,m3     Bottom row of the Wigner3j symbol; positive orders
% method       'automatic' (default), or force:
%              'gl'    By Gauss-Legendre integration
%              'recu'  By recursion formula
%
% OUTPUT:
%
% W            The result, perhaps with all the recursions as a vector.
%              This makes it hard to save and retrieve. As this was
%              changed halfway, if something goes wrong, recalculate the
%              saved ones! 
%
% EXAMPLE:
%
% wigner3j('demo1')
% l=round(100*rand(1)); [l wigner3j(l,l,0,0,0,0)^2*(2*l+1)]
%
% Last modified by fjsimons-at-alum.mit.edu, 28.3.2005

defval('m1',0)
defval('m2',0)
defval('m3',0)
defval('method','automatic')

if (m1<0 | m2<0 | m3<0) 
  % Only special cases are allowed, with the Condon-Shortley phase
  % Top row is of the form (l+a b l) with a=0,1,2, and b=1,2
  % Bottom row is of the form (-m 0 m)
  if l2==1 & m2==0 & m1==-m3 & l1==l3+1
    l=l3; m=abs(m1);
    W=(-1)^(l+m+1)*sqrt(2*((l+1)^2-m^2)/((2*l+3)*(2*l+2)*(2*l+1)));
    return
  end
  if l2==2 & m2==0 & m1==-m3 & l1==l3+2
    l=l3; m=abs(m1);
    W=(-1)^(l+m)*sqrt(6*((l+1)^2-m^2)*((l+2)^2-m^2)/...
                        ((2*l+5)*(2*l+4)*(2*l+3)*(2*l+2)*(2*l+1)));
    return
  end
  if l2==2 & m2==0 & m1==-m3 & l1==l3+1
    l=l3; m=abs(m1);
    W=(-1)^(l+m+1)*2*m*sqrt(6*((l+1)^2-m^2)/...
			    ((2*l+4)*(2*l+3)*(2*l+2)*(2*l+1)*2*l));
    return
  end
  if l2==2 & m2==0 & m1==-m3 & l1==l3
    l=l3; m=abs(m1);
    W=(-1)^(l+m)*2*(3*m^2-l*(l+1))/...
                sqrt((2*l-1)*(2*l)*(2*l+1)*(2*l+2)*(2*l+3));
    return
  else
    error('This don''t work so good.')
  end
end

if ~isstr(l1)
  fnpl=sprintf('%s/W3J-%i-%i-%i-%i-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'WIGNER'),l1,l2,l3,m1,m2,m3);
  if exist(fnpl,'file')==2 % & 1==3
    eval(sprintf('load %s',fnpl))
%    disp(sprintf('%s loaded by WIGNER3J',fnpl))
  else
    % Negative m's are allowed unless the result is computed by integration 
    if abs(m1)>l1 | abs(m2)>l2 | abs(m3)>l3
      error('Mismatch of angular degrees and orders')
    end 

    % Conditions of Dahlen and Tromp Eq. (C.186) % BUT NO NEGATIVE M
    if ~triangle(l1,l2,l3) | m1~=(m2+m3)
      W=0;
      return
    end

    % Conditions of Dahlen and Tromp Eq. (C.219)
    if m1==0 & m2==0 & m3==0 & mod(l1+l2+l3,2)
      W=0;
      return
    end

    % Via Gauss-Legendre integration, forced or by choice
    if (l1~=l2 | m1~=0 | m2~=0 | m3~=0 | strcmp(method,'gl') )
      % If the m is odd then the associated is not a polynomial
      % and this is not the best way; would need to add new recursion
      % relations to cover this case.
      integrand=inline(sprintf(...
	  ['rindeks(legendre(%i,x),%i).*',...
	   'rindeks(legendre(%i,x),%i).*',...
	   'rindeks(legendre(%i,x),%i)'],...
	  l1,m1+1,l2,m2+1,l3,m3+1));
      % Sign information is lost, however.
      % Dahlen and Tromp Eq. (C.226)
      W=sqrt(gausslegendre([-1 1],integrand,l1+l2+l3)/2);
      %disp('Wigner3j by Gauss-Legendre integration')
      eval(sprintf('save %s W',fnpl))  
      return
    end

    % Dahlen and Tromp Eqs. (B.49) and (C.226)
    if l3==0 & m1==m2 & ~strcmp(method,'recu')
      % The triangle equality with l3=0 implies l1==l2
      % The initial check on m implies that m3==0
      if (abs(l1)+abs(m1))<21
	%disp('Wigner3j by normalization formula Dahlen and Tromp B.49')
	W=sqrt(factorial(l1+m1)/factorial(l1-m1)/(2*l1+1));
	eval(sprintf('save %s W',fnpl))  
	return
      else
	% If factorial fails need to use Gauss-Legendre
	W=wigner3j(l1,l2,l3,m1,m2,m3,'gl');
	return
      end
    end

    % All-zero bottom row; Dahlen and Tromp Eq. (C.220); sign positive
    if m1==0 & m2==0 & m3==0 & (l1+l2+l3)<20 & ~strcmp(method,'recu')
      S=(l1+l2+l3); 
      %disp('Wigner3j by summation formula Dahlen and Tromp C.220')
      W=sqrt(1/factorial(S+1)*...
	     factorial(S-2*l1)*factorial(S-2*l2)*factorial(S-2*l3))*...
	factorial(S/2)/...
	factorial(S/2-l1)/factorial(S/2-l2)/factorial(S/2-l3);
      eval(sprintf('save %s W',fnpl))  
      return
    end

    % All-zero bottom row and l1==l2; 
    % Recursion http://functions.wolfram.com/07.39.17.0013.01
    if m1==0 & m2==0 & m3==0 & l1==l2
      % By upward recursion and saving all the terms
      W=repmat(NaN,1,l3+1);
      % Dahlen and Tromp Eq. (B.49)
      W(1)=1/sqrt(2*l1+1); % Zeroth term
      for j=2:2:l3 
	% All odd terms zero from triangle rule with l1==l2 
	W(j)=0;
	% Keep sign positive
	W(j+1)=(j-1)/j*...
	       sqrt((4*l1^2+4*l1-j^2+2*j)/(4*l1^2+4*l1-j^2+1))...
	       *W(j-1);
      end
      %disp('Wigner3j by recursion formula Mathematica')
      eval(sprintf('save %s W',fnpl))  
    end
  end  
elseif strcmp(l1,'demo1')
  j=0; m=0;
  for l=0:100
    d(l+1)=wigner3j(l,l,j,m,m,m)^2-legendreprodint(l,m,l,m,-1,'gl')/2;
  end 
  clf
  plot(0:100,abs(d)/eps,'LineW',2)  
  xlabel('l')
  ylabel(sprintf('W_{ll}^0 difference %s eps','\times'))
  title('Difference between WIGNER3J and direct GL integration');
  grid on
elseif strcmp(l1,'demo2')
  m=0; l=round(rand*100);
  for j=0:100
    D=wigner3j(l,l,j,m,m,m);
    d(j+1)=D(end)-wigner3j(l,l,j,m,m,m,'gl');
    clear D
  end
  clf
  plot(0:100,abs(d)/eps,'LineW',2)  
  xlabel('j')
  ylabel(sprintf('W_{%i%i}^j difference %s eps',l,l,'\times'))
  title('Difference between WIGNER3J using recursion and GL integration');
  grid on
else
  error('Specify valid option')
end
