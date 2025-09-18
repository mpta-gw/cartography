function y = normCalPriorComplex(Lambda, phi, omega, omega_opt, ...
                                 sigma, delta_Lambda, delta_phi)
%  y = normCalPriorComplex(Lambda, phi, omega, omega_opt, ...
%                          sigma, delta_Lambda, delta_phi)
%
%  is the integrand for numerical integration in
%  omegaPosteriorComplex().  omega_opt and sigma are the measured
%  (complex) Omega estimate and corresponding errorbar; omega is
%  a hypothetical actual value for Omega; Lambda is a vector of
%  log-amplitude values and phi a single phase value for the
%  calibration.  delta_Lambda and delta_phi are the standard
%  deviations of the corresponding gaussians.
%
%  "help omegaPosteriorComplex" for more details.
%
%  Routine written by John T. Whelan
%  Contact john.whelan@ligo.org
%
% $Id: normCalPriorComplex.m,v 1.1 2006-04-14 23:45:30 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = exp(-Lambda.^2/(2*delta_Lambda^2)) ...
    * exp(-phi^2/(2*delta_phi^2)) ...
    .* exp( - abs( omega_opt - omega*exp(-Lambda-1i*phi) ).^2 ...
            / (2*sigma^2) );
return
