function [omega_posterior, omega, likelihood] ...
    = omegaPosteriorComplex(ptestimate, sigma, ...
                            omega_min, omega_max, omega_prior, ...
                            doMarginalize, ...
                            amp_uncertainty_pct, phase_uncertainty_deg, ...
                            Nsigma_Lambda, Nsigma_phi, ...
                            intMethod, num_Lambda, num_phi);
%
%  [omega_posterior, omega, likelihood] ...
%     = omegaPosteriorComplex(ptestimate, sigma, ...
%                             omega_min, omega_max, omega_prior, ...
%                             doMarginalize, ...
%                             amp_uncertainty_pct, phase_uncertainty_deg, ...
%                             Nsigma_Lambda, Nsigma_phi, ...
%                             intMethod, num_Lambda, num_phi);
% 
%  calculates a posterior PDF from the specified complex point estimate
%  and errorbar, optionally marginalizing over the calibration
%  uncertainty, which is taken to be gaussian-distributed in
%  log-amplitude and phase.
%  
%  Inputs:
%
%    ptestimate = the complex point estimate of Omega
%
%    sigma = the one-sigma error bar associated with the real
%               and imaginary parts of the point estimate
%
%    omega_min = the lower end of the range of Omegas
%
%    omega_max = the upper end of the range of Omegas
%
%    num_omega = prior on Omega as a vector of relative probabilities for
%                evenly spaced Omega values; need not be normalized
%
%    doMarginalize = boolean flag indicating whether or not to
%                    marginalize the likelihood function over
%                    calibration uncertainty.  (Subsequent
%                    arguments are irrelevant and unneeded if
%                    doMarginalize = false.)  Default = false
%
%    amp_uncertainty_pct = percent uncertainty in the calibration
%                          amplitude, expressed as a one-sigma
%                          error on log-amplitude
%
%    phase_uncertainty_deg = one-sigma uncertainty in the
%                            calibration phase in degrees,
%
%    Nsigma_Lambda = number of standard deviations in log-amplitude
%                    below and above zero for range of numerical
%                    integral. Default = 3. 
%
%    Nsigma_phi = number of standard deviations in phase below and
%                 above zero for range of numerical integral.
%                 Default = 3.
%
%    intMethod = integration method to use for numerical integral
%                over calibration amplitude and phase.  Alowed
%                values are 'brute' (brute-force sum over
%                fixed-resolution grid; slow but reliable) and
%                'quad' (matlab dblquad function; faster but less
%                robust).  (Subsequent arguments are irrelevant and
%                unneeded unless intMethod = 'brute'.)
%                Default = 'quad'.
%
%    num_Lambda = number of points in log-amplitude direction for
%                 brute-force numerical integration.  Default = 100.
%
%    num_phi = number of points in phase direction for brute-force
%              numerical integration.  Default = 100.
%
%  Outputs:
%
%    omega_posterior = posterior PDF on Omega, normalized so that its
%              integral from omega_min to omega_max is unity
%
%    omega = a vector of Omega values corresponding to the vector
%            omega_posterior
%
%    likelihood = the unnormalized likelihood function
%
%  This routine calculates a likelihood function for the actual
%  observed real and imaginary cross-correlation estimates of Omega
%  and the associated error bar, for a set of hypothetical actual
%  values of Omega.
%
%  If doMarginalize = false, the unnormalized likelihood function
%  is just a Gaussian
%  $$
%  P(xy|\Omega\sigma) \propto \exp(-[x-\Omega]^2/2\sigma^2)
%  $$
%  where x+iy is the measured Omega, \sigma is the corresponding
%  errorbar, and \Omega is the unknown true value of Omega.  (The factor
%  exp(-y^2/2\sigma^2) does not contain \Omega and would thus drop out of the
%  normalized prior, so it is omitted.
%
%  If doMarginalize = true, the unnormalized likelihood function
%  is constructed by a numerical integral over the log-amplitude
%  and phase of the calibration error
%  $$
%  P(xy|\Omega\sigma) \propto \iint d\Lambda d\phi
%  \exp( - |x+iy-\Omega\exp[\Lambda+i\phi]|^2/2\sigma^2
%        - \Lambda^2/2\sigma_\Lambda^2 - \phi^2/2\sigma_\phi^2 )
%  $$
%  The integrand is calculated with the function normCalPriorComplex()
%
%  The likelihood function is then restricted to the range from
%  omega_min to omega_max, multiplied by the specified prior, and
%  normalized so its integral (NOT sum) over that range is unity.
%
%  Routine written by John T. Whelan
%  (derived from bayesianMarginalize by Joseph D. Romano)
%  Contact john.whelan@ligo.org
%
% $Id: omegaPosteriorComplex.m,v 1.3 2006-04-16 17:16:15 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default params
try
  doMarginalize;
catch
  doMarginalize = false;
end;
try
  Nsigma_Lambda;
catch
  Nsigma_Lambda = 3;
end;
try
  Nsigma_phi;
catch
  Nsigma_phi = 3;
end;
try
  intMethod;
catch
  intMethod = 'quad';
end;
try
  num_Lambda;
catch
  num_Lambda = 100;
end;
try
  num_phi;
catch
  num_phi = 100;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure omega_prior is a column vector
omega_prior = omega_prior(:);

% discrete omega values
num_omega = length(omega_prior);
omega = linspace(omega_min, omega_max, num_omega);
omega = transpose(omega);
delta_omega = omega(2)-omega(1);
                                                                                
% initialize arrays
likelihood = zeros(num_omega,1);
omega_posterior = zeros(num_omega,1);

if doMarginalize

  % calib range 
  delta_Lambda = log(1+amp_uncertainty_pct/100);
  Lambda_max = Nsigma_Lambda * delta_Lambda;
  Lambda_min = -Lambda_max;
  delta_phi = phase_uncertainty_deg * pi / 180;
  phi_max = min(Nsigma_phi * delta_phi, pi);
  phi_max = delta_phi;
  phi_min = -phi_max;

  % discrete Lambda and phi for integral
  Lambda = linspace(Lambda_min, Lambda_max, num_Lambda);
  Lambda = transpose(Lambda);
  dLambda = Lambda(2)-Lambda(1);
  phi = linspace(phi_min, phi_max, num_phi);
  phi = transpose(phi);
  dphi = phi(2)-phi(1);

  switch lower(intMethod)
   case 'brute';
    % loop over omega values
    for ii=1:num_omega
      if (mod(ii,100) == 0)
        fprintf('%d out of %d integrations finished\n', ii, num_omega);
      end;
      integral = 0;
      for jj=1:num_phi
        int = normCalPriorComplex(Lambda, phi(jj), omega(ii), ptestimate, ...
                                  sigma, delta_Lambda, delta_phi);      
        integral  = integral + sum(int) * dLambda * dphi;
      end;
      
      % marginalised likelihood (not normalised however)
      likelihood(ii) = integral;
    end; % for ii=1:num_omega
   case 'quad'
    % loop over omega values
    for ii=1:num_omega
       
%       if (mod(ii,100) == 0)
%         fprintf('%d out of %d integrations finished\n', ii, num_omega);
%       end;
      % integrate over calibration values
      integral = dblquad(@normCalPriorComplex, ...
                         Lambda_min, Lambda_max, phi_min, phi_max, ...
                         [], [], omega(ii), ptestimate, sigma, ...
                         delta_Lambda, delta_phi);
      
      % marginalised likelihood (not normalised however)
      likelihood(ii) = integral;
    end; % for ii=1:num_omega
   otherwise
    error('Unrecognized integration method');
  end; % switch lower(intMethod)
  else
    % non-marginalised likelihood (not normalised however)
    likelihood = exp(-0.5*(((real(ptestimate)-omega)/sigma ).^2));
end

% calculate posterior distribution (truncated to omega_min to omega_max)
omega_posterior = likelihood .* omega_prior;

% calculate normalisation factor
norm = 1/(delta_omega*sum(omega_posterior));

omega_posterior = norm*omega_posterior;

return
