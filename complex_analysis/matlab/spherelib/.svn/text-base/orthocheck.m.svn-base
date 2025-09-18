function [ngl1,ngl2,com,Vc]=orthocheck(C,V,TH,m,N)
% [ngl1,ngl2,com,Vc]=ORTHOCHECK(C,V,TH,m,N)
%
% Checks the orthonormality of a SINGLE-ORDER spherical harmonic
% expansion over the UNIT SPHERE and over a SINGLE or DOUBLE spherical CAP.
% If no EIGENVALUES are known, calculates them by integration.
% All in fully normalized real harmonics (to 1 over unit sphere).
%
% INPUT:
%
% C        Spherical harmonic expansion coefficients
% V        Eigenvalues, in a matrix, if known
% TH       Radius of the spherical polar cap, in RADIANS
% m        Angular order, m>=0
% N        1 Single polar cap [default]
%          2 Double polar cap
%
% OUPUT:
%
% ngl1     Number of GL points on the unit sphere to satisy orthogonality
% ngl2     Number of GL points on the domain to satisy orthogonality
% com      Center of mass of SINGLE cap (see Freeden/Narcowich)
% Vc       Eigenvalues deduced from the integration
%
% Displays a message with the mean absolute error over the coefficients,
% the unit sphere, and the domain. The number of GL points are obtained
% by assuming that the orthogonality is satisfied and then determining the
% number which best verifies this, which is useful for non-zero order.
%
% Last modified by fjsimons-at-alum.mit.edu, 16.08.2005

defval('N',1)

tol=10^-12;

L=size(C,1)+abs(m)-1;

% Verify normalization of the coefficients; this should be zero
% This here is testing SDW (2005) equation (4.8), first part
% This here is testing SDW (2005) equation (5.9), first part
err(1)=mean(mean(abs(C'*C-eye(size(C,2)))));
% Hmm... why not also test C'*K*C'? Should be the eigenvalues.

% Verify normalization to unity over the unit sphere
% This here is testing SDW (2005) equation (4.9), first part
% Upping the number of GL points is crucial
if m==0
  % At least two integration points, to cover the L=0 case without messing
  % with pl(m)2th
  ngl=max(2*L,2);
  % Normalization over the unit sphere
  [w,x]=gausslegendrecof(ngl,[],[-1 1]);
  % This is SDW (2005) equation (5.10) at the GL points
  % BUT IT ALSO INCLUDES sqrt(2-dom) of equation (5.12)
  EGL=pl2th(C,acos(x),1);
  % This should be zero
  err(2)=orvo1(w,EGL,m);
  ngl1=ngl;
  if N==1
    % Center of mass
    com=2*pi*diag(EGL'*diag(w.*x(:))*EGL);
    com=com';
  else
    com=0;
  end
elseif m~=0
  ngl=2*L+[1:2:100];
  % Iteratively determines the best number of integration points
  for index=1:length(ngl)
    [w,x]=gausslegendrecof(ngl(index),[],[-1 1]);
    % This is SDW (2005) equation (5.10) at the GL points
    % BUT IT ALSO INCLUDES sqrt(2-dom) of equation (5.12)
    EGL=plm2th(C,acos(x),m,1);
    trerr(index)=orvo1(w,EGL,m);
    if trerr(index)<tol; break; end
  end
  [err(2),jk]=min(trerr);
  ngl1=ngl(index);
  com=NaN;
  if index==length(ngl); 
    disp('Could do with more integration accuracy')
  end
end

clear trerr

% Verify normalization to eigenvalue over the polar cap(s)
% This here is testing SDW (2005) equation (4.9), second part
% If ngl is too big leads to inaccurate results!
% If V is a matrix, it's the standard case, and we check.
% If V is a vector then it's the Sturm-Liouville eigenvalue and we
% need to still calculate the localization eigenvalue.
intv=[cos(TH) 1];
if m==0
  [w,x]=gausslegendrecof(ngl,[],intv);
  % This is SDW (2005) equation (5.7) at the GL points
  % BUT IT DOES INCLUDE sqrt(2-dom) of equation (5.9)
  EGL=pl2th(C,acos(x),1);
  if N==1 % SIMPLE CAP
    [err(3),Vc]=orvo2(w,EGL,m,V);
  elseif N==2 % DOUBLE CAP
    % The top half
    % The sqrt(2) had already been taken care of, thus the longitudinal
    % integral is either pi or 2pi
    orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
    % Add the bottom half, same weight, flipped x
    EGL=pl2th(C,acos(-flipud(x)),1);
    % The sqrt(2) had already been taken care of, thus the longitudinal
    % integral is either pi or 2pi
    orv4=orv4+[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
    if prod(size(V))~=length(V) & ~isempty(V)
      % Then, not from GRUNBAUM2, eigenvalues known
      err(3)=mean(mean(abs(orv4-V)));
      Vc=V;
    else
      % Add it in for GRUNBAUM2
      Vc=diag(orv4)';
      % But the off-diagonal terms must be zero for orthogonality
      err(3)=mean(mean(abs(orv4-diag(Vc))));
    end
  end
  ngl2=ngl;
elseif m~=0
  % Iteratively determines the best number of integration points
  for index=1:length(ngl)
    [w,x]=gausslegendrecof(ngl(index),[],intv);
  % This is SDW (2005) equation (5.7) at the GL points
  % BUT IT DOES INCLUDE sqrt(2-dom) of equation (5.9)  
  EGL=plm2th(C,acos(x),m,1);
    if N==1 % SIMPLE CAP
      [trerr(index),Vc]=orvo2(w,EGL,m,V);
      if trerr(index)<tol; break; end
    elseif N==2 % DOUBLE CAP
      % The top half
      % The sqrt(2) had already been taken care of, thus the longitudinal
      % integral is either pi or 2pi
      orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
      % Add the bottom half, same weight, flipped x
      EGL=plm2th(C,acos(-flipud(x)),m,1);
      % The sqrt(2) had already been taken care of, thus the longitudinal
      % integral is either pi or 2pi
      orv4=orv4+[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
      if prod(size(V))~=length(V) & ~isempty(V)
	% Then, not from GRUNBAUM2, eigenvalues known
	trerr(index)=mean(mean(abs(orv4-V)));
	Vc=V;
	if trerr(index)<tol; break; end
      else
	% Add it in for GRUNBAUM2
	Vc=diag(orv4)';
	% But the off-diagonal terms must be zero for orthogonality
	trerr(index)=mean(mean(abs(orv4-diag(Vc))));
	if trerr(index)<tol; break; end
      end
    end
  end
  [err(3),jk]=min(trerr);
  ngl2=ngl(index);
  if index==length(ngl); 
    disp('Could do with more integration accuracy')
  end
end

if any(err>tol)
  warning(sprintf('Normalization criteria NOT satisfied; mean errror %8.3e',...
	       mean(err)))
else
  disp(sprintf('Normalization criteria satisfied to %8.3e',...
	       mean(err)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function erro=orvo1(w,EGL,m)
% Calculate the latitudinal and longitudinal integral
% For the diagonal elements only?
% orv1=w(:)'*[EGL.^2]*(1+(m==0))*pi;
% For the entire matrix
% The sqrt(2) had already been taken care of, thus the longitudinal
% integral is either pi or 2pi
orv2=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
erro=mean(mean(abs(orv2-eye(size(orv2)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [erro,Vc]=orvo2(w,EGL,m,V)
% Calculate the latitudinal and longitudinal integral
% For the diagonal elements only?
% orv3=w(:)'*[EGL.^2]*(1+(m==0))*pi;
% For the entire matrix
% The sqrt(2) had already been taken care of, thus the longitudinal
% integral is either pi or 2pi
orv4=[EGL'*diag(w)*EGL]*(1+(m==0))*pi;
if prod(size(V))~=length(V) & ~isempty(V)
  erro=mean(mean(abs(orv4-V)));
  Vc=V;
else
  % We don't check but calculate the eigenvalues as the diagonal
  Vc=diag(orv4)';
  % But the off-diagonal terms must be zero for orthogonality
  erro=mean(mean(abs(orv4-diag(Vc))));
end

