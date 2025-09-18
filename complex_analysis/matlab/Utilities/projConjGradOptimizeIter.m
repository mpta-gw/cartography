function [x0,CG]=projConjGradOptimizeIter(x0,CG)

% calculates one conjugate gradient iteration towards optimizing
%  y = (func)(x)
% While x can be n-dim, y has to be 
%
% uses
% NUMTOOLS.stepsize
% NUMTOOLS.func
% NUMTOOLS.lastx0
% NUMTOOLS.lasty0

global NUMTOOLS

try, CG; catch try, CG=NUMTOOLS.CG; catch CG=[]; end; end;
try, x0; catch x0=NUMTOOLS.lastx0; end;

  if length(CG)==0
    [y,G]=calcGradientJacobi(x0);
    if length(y)>1
      % (func) has to be 1-dim - use norm instead
      NUMTOOLS.alternateFunc = NUMTOOLS.func;
      NUMTOOLS.func          = 'normOfFunction';
      [y,G]=calcGradientJacobi(x0);
    end
    CG.g  = -transpose(G);
    CG.xi = CG.g;
  end  



  [f,fp,fpp]=calcdiffdiff(x0,CG.xi);
  dx0=-(fp)/(fpp)*CG.xi/norm(CG.xi);
  x0=x0+dx0*.5;
  [y,G]=calcGradientJacobi(x0);

  h      = CG.xi;
  CG.xi  = transpose(G);
  gam    = (transpose(CG.xi  )*CG.xi) / (transpose(CG.g)*CG.g);
  CG.g   = -CG.xi;
  CG.xi  = CG.g + gam.*h;
      
  CG.err          = norm(G);
  NUMTOOLS.lastx0 = x0;
  NUMTOOLS.lasty0 = y;
  NUMTOOLS.CG     = CG;
