function [mode, mu_UL, bayesfactor] = ...
  bayesianUL(xbar, sigma, confidence, mu_min, mu_max)
%
%  bayesianUL --- calculates mode, Bayesian UL, and Bayes factor
%
%  bayesianUL(xbar, sigma, confidence, mu_min, mu_max) calculates
%  the mode, Bayesian UL, and Bayes factor for a Gaussian distributed 
%  pdf (mean xbar, stddev sigma--assumed known) truncated between 
%  mu_min and mu_max.  The bayes factor is the ratio of the
%  probability of the data assuming no signal to the probability of 
%  the data assuming a signal model described by a single parameter.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%  
%  $Id: bayesianUL.m,v 1.8 2005-10-11 14:34:34 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default value for mu_min=0
try 
  mu_min;
catch
  mu_min=0;
end

% default value for mu_max=inf
try 
  mu_max;
catch
  mu_max=inf;
end

% calculate mode 
mode = xbar;
mode(find(xbar < mu_min)) = mu_min;

% calculate upper limit
mu_UL = xbar + sqrt(2)*sigma.*...
                erfcinv( (1-confidence)*erfc((mu_min-xbar)./(sqrt(2)*sigma))...
                          + confidence *erfc((mu_max-xbar)./(sqrt(2)*sigma)) );

% calculate bayes factor
S12 = (1/sqrt(2*pi))*(1./sigma)*2.*(mu_max-mu_min)./...
              ( erfc((mu_min-xbar)./(sqrt(2)*sigma)) - ...
                erfc((mu_max-xbar)./(sqrt(2)*sigma)) );
R12 = exp(-0.5*(xbar./sigma).^2);
bayesfactor = S12.*R12;

return
