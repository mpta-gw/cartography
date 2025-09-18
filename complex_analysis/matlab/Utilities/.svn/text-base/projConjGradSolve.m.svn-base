function x0=projConjGradSolve(x0,maxerr,maxiter,stepsize,func);

% iteratively calls [x0,CG]=conjGradOptimizeIter(x0,CG) to optimize
%  (func)(x0)
% Note that the output of func has to be 1-dimensional
%
% uses
% NUMTOOLS.stepsize
% NUMTOOLS.func
% NUMTOOLS.maxerr
% NUMTOOLS.maxiter
% NUMTOOLS.lastx0
% NUMTOOLS.lasty0

global NUMTOOLS

try, NUMTOOLS.func     = func;     end;
try, NUMTOOLS.stepsize = stepsize; end;
try, NUMTOOLS.maxiter  = maxiter;  end;
try, NUMTOOLS.maxerr   = maxerr;   end;

try, NUMTOOLS.maxiter; catch NUMTOOLS.maxiter=1000; end;
try, NUMTOOLS.maxerr;  catch NUMTOOLS.maxerr =1e-6; end;
 
try, x0; catch x0=zeros(size(NUMTOOLS.stepsize)); end;

ii=0;
err=NUMTOOLS.maxerr;
CG=[];
NUMTOOLS.VALUES=[];
while  ii<NUMTOOLS.maxiter && err>=NUMTOOLS.maxerr,
  [x0,CG]=projConjGradOptimizeIter(x0,CG);
  x0=projectPositive(x0);
  val=minimizeThis(x0)-NUMTOOLS.BESTVALUE;
  NUMTOOLS.VALUES(ii+1)=val;
  display(['******************* Offset from absolute minimum ',num2str(val)]);
  err=CG.err;
  ii=ii+1;
  disp(['Iteration ',num2str(ii),'; err = ',num2str(err)]);
end


