function y = convertResponse(fx, x, flow, deltaF, numFreqs, invert, powerOfTen)
%
%  convertResponse --- interpolates, inverts, or scales a response function
%
%  convertResponse(fx, x, flow, deltaF, numFreqs, invert, powerOfTen)
%  - interpolates a response function x originally evaluated at the 
%    discrete frequencies fx to a new set of discrete frequencies:
%    f = flow + deltaF*[0:numFreqs-1]
%  - inverts from e.g., strain/count to count/strain (if invert==1)
%  - scales by 10^powerOfTen to get e.g., count/attostrain 
%
%  The result is output to a frequency-series structure.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: convertResponse.m,v 1.2 2005-02-24 13:48:17 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct desired frequencies
f = flow + deltaF*[0:numFreqs-1]';

% interpolate (real and imag parts separately)
x_real = real(x);
x_imag = imag(x);

y_real = interp1(fx, x_real, f, 'linear');
y_imag = interp1(fx, x_imag, f, 'linear');
y_temp = y_real + sqrt(-1)*y_imag;

% invert response (if desired)
if invert == 1
  y_temp = 1./y_temp;
end

% scale response
y_temp = (10^powerOfTen)*y_temp;

% fill structure
y.data = y_temp;
y.flow = flow;
y.deltaF = deltaF;

return

