function [response, responseOK] = ...
  calculateResponse(t, f, R0, C0, alpha, gamma, epoch, channeltype)
%
%  calculateResponse --- calculates the instrument response function 
%
%  calculateResponse(t, f, R0, C0, alpha, gamma, epoch, channeltype)
%  calculates the values of the (time-dependent) instrument response
%  function at a particular epoch from the (time-dependent) calibration
%  information.
%  
%  R0 and C0 are arrays containing the values of the reference response 
%  function and sensing function (cavity gain) at t=0 (evaluated at 
%  discrete frequencies f)
%
%  alpha and gamma are arrays containing the values of the cavity factor 
%  (C = alpha*C0) and open loop factor (H = gamma*[C0*R0 - 1] ) evaluated
%  at discrete times t.
%  
%  epoch is the GPS time of the calculated response function.
%
%  channeltype is either 'AS_Q' or 'DARM_ERR' (if omitted, 'AS_Q' is assumed)
%  The 'LSC-' is optional.
%
%  The output is an array containing the values of the response function
%  R(f,epoch) evaluated at the discrete frequencies f.
%
%  if channeltype is 'AS_Q' the response function is
%
%   R(f,epoch) = [ 1 + gamma(epoch)*(C0(f)*R0(f) - 1) ] 
%                --------------------------------------
%                      alpha(epoch)*C0(f)
%
%  if channeltype is 'DARM_ERR' the response function is
%
%   R(f,epoch) = [ 1 + gamma(epoch)*(C0(f)*R0(f) - 1) ] 
%                --------------------------------------
%                      gamma(epoch)*C0(f)
%
%
%  responseOK is either false or true depending on whether or not 
%  calibration data existed at the desired epoch
% 
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: calculateResponse.m,v 1.7 2005-05-11 15:51:28 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that frequency-series have the same length
if ( length(R0) ~= length(C0) )
  error('size mismatch');
end

% check that times-series have the same length
if ( length(alpha) ~= length(gamma) )
  error('size mismatch');
end

% check for bad epoch
if ( epoch < t(1) )
  error('desired epoch before calibration line data was taken');
end;
if ( epoch > t(end) )
  error('desired epoch after calibration line data was taken');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find calibration time closest to epoch
[ignore, i] = min(abs(t-epoch));

% set alpha and gamma values appropriately
t_epoch = t(i);
alpha_epoch = alpha(i);
gamma_epoch = gamma(i);

if (alpha_epoch < 1e-6)

  fprintf('no calibration data available for epoch %d\n', t_epoch);
  response = 0;
  responseOK = false;

else

  % if DC component of sensing function = 0, set C0(DC)=C0(deltaF)
  if ( (f(1) == 0) & (C0(1) == 0) )
    fprintf('NOTE!!! DC value of sensing function = 0\n');
    fprintf('Setting DC component to its value at deltaF in order to calculate response function\n');
    C0(1)=C0(2);
  end

  % calculate response
  try
    channeltype;
  catch
    channeltype = 'AS_Q';
  end;
  if strncmp(channeltype, 'LSC-', 4)
    channeltype = channeltype(5:end);
  end;

  if strncmp(channeltype, 'AS_Q', length(channeltype))
    response = (1+gamma_epoch*(C0.*R0 - 1))./(alpha_epoch*C0);
    responseOK = true;
  elseif strncmp(channeltype, 'DARM_ERR', length(channeltype))
    response = (1+gamma_epoch*(C0.*R0 - 1))./(gamma_epoch*C0);
    responseOK = true;
  elseif strncmp(channeltype, 'ETMX_EXC_DAQ', length(channeltype))
    response = R0;
    responseOK = true;
  else
    errmsg = sprintf('Unrecognized channel type %s',channeltype);
    error(errmsg);
  end;
end

return

