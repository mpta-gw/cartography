function R=timeRotation(Lmax,siderealTime,diagonly)

% performs a rotation around the earth axis by siderealTime
% Example: To get the gamma_lm_t coefficients at a given time:
%
%            glm = calGammaLM(gammaLM_coeffsPath,detPair, lmax, numFreqs, flow, deltaF);
%            gamma_lm_t=timeRotation(lmax,siderealTime)*glm.data;
%
% If length(siderealTime)=1, this function is equivalent to
%    zRotation(Lmax,-pi/12*siderealTime,diagonly)
%
% Since R is a diagonal Matrix, there is an option to just return the diagonal.
% Default is: R is rotation matrix for length(siderealTime)=1
%             R is #lm x #timeBins matrix if length(siderealTime)>1

if length(siderealTime)>1,
  diagonly=true;
else
  try, diagonly; catch diagonly=false; end;
end;

[lvec,mvec]=getLMvec(Lmax);

% Sign: rotate coefficients in opposite direction
r=exp(i*pi*mvec*siderealTime/12);

if diagonly
  R=r;
else
  R=diag(r);
end;
