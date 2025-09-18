function [s] = KSStatistic(d)
%
%  KSStatistic(d) returns the value of the Kolmogorov-Smirnov 
%  statistic for input d.  (Refer to Eq. 14.3.7, p.624 in Numerical 
%  Recipes in C 2nd Ed.)  To prevent underflow, the infinite sum is 
%  truncated at a value k for the summation index for which 
%  exp(-2*(k*d)^2) = exp(-10).
%
%  Routine written by Albert Lazzarini.
%  Contact lazz@ligo.caltech.edu
%
%  $Id: KSStatistic.m,v 1.1 2005-10-19 17:51:01 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imax = floor(10.0/(sqrt(2)*d) + 1);

s = 0.;
for i = 1:imax,
    s = s + 2*(-1)^(i + 1)*exp(-2*i^2*d^2);
end

return
