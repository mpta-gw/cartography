function varargout=calcdiffdiff(x0,dx)

% calculates the 1st and 2nd derivative of (func) at x0 into the direction dx
%
% returns:
%     y     : Value at x0
%     y1    : 1st derivative into direction dx
%     y2    : 2nd derivative into direction dx
%
% uses
% NUMTOOLS.stepsize
% NUMTOOLS.func

global NUMTOOLS

dims=length(x0);

y=eval([NUMTOOLS.func,'(x0)']);
varargout{1}=y;

if nargout>1  
  % How far should I go along dx?
  dx=dx./ max(max(abs(dx./NUMTOOLS.stepsize)));
  adx=norm(dx);
  
  yp=eval([NUMTOOLS.func,'(x0+dx)']);
  ym=eval([NUMTOOLS.func,'(x0-dx)']);
  y1=(yp-ym)/(2*adx);
  varargout{2}=y1;
  if nargout>2
    y2=(yp-2*y+ym)/(adx^2);
    varargout{3}=y2;
  end
end

return



