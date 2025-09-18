function varargout=calcGradientJacobi(x0)

% calculates the gradient of (func) at x0
% (x0 has the size Mx1)
%
% returns:
%     y    : Value of (func) at x0       (size Nx1)
%     G    : Gradient of (func) at x0    (size NxM)
%     J    : Gradient of (func) at x0    (size MxMxN),
%            i.e. if y is more than 1-dim, the third dimension is used
%
% uses
% NUMTOOLS.stepsize
% NUMTOOLS.func

global NUMTOOLS

dims=length(x0);

y=eval([NUMTOOLS.func,'(x0)']);

if nargout>1
  dimsy=length(y);
  G=zeros(dimsy,dims);
  if nargout>2, J=zeros(dims,dims,dimsy); end;
  for ii=1:dims
    dx=zeros(dims,1);
    dx(ii)=NUMTOOLS.stepsize(ii);
    yp=eval([NUMTOOLS.func,'(x0+dx)']);
    ym=eval([NUMTOOLS.func,'(x0-dx)']);
    G(:,ii)=(yp-ym)/(2*dx(ii));
    if nargout>2
      J(ii,ii,:)=(yp-2*y+ym)/(dx(ii).^2);
      for jj=1:(ii-1)
        dx2=zeros(dims,1);
   	dx2(jj)=NUMTOOLS.stepsize(jj);
	ypp=eval([NUMTOOLS.func,'(x0+dx+dx2)']);
	ymp=eval([NUMTOOLS.func,'(x0-dx+dx2)']);
	ypm=eval([NUMTOOLS.func,'(x0+dx-dx2)']);
	ymm=eval([NUMTOOLS.func,'(x0-dx-dx2)']);
        J(ii,jj,:)=(ypp-ymp-ypm+ymm)/4/dx(ii)/dx2(jj);
	J(jj,ii,:)=J(ii,jj,:);
      end
    end
  end
end

switch nargout
  case 1,
      varargout{1}=y;
  case 2,
      varargout{1}=y;
      varargout{2}=G;
  case 3,
      varargout{1}=y;
      varargout{2}=G;
      varargout{3}=J;
end
return



