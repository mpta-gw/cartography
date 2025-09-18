function CP = SpH2CrossPower(glm, plm, siderealTime, numFreqs, flow, deltaF, H, beta)

% calculates the cross-power due to a background given by H*plm
% input: glm          struct from calGammaLM
%        plm          spherical harmonics coefficients for the P(Omega)
%        siderealTime evaluate at this sidereal time
%        numFreqs     number of frequency bins
%        flow         lowest frequency
%        deltaF       frequency bin separation
%        H            power spectrum, either vector with numFreq elements,
%                     or value at 100Hz (scalar)
%        beta         power law (if H is scalar), default =0
%
% output:CP           Cross-Power

try, beta; catch beta=0; end
freq = glm.flow + (0:glm.numFreqs-1)'*glm.deltaF;
if length(H)==1
  H=H*((freq/100.0) .^ beta);
end

%CP=  ( transpose(glm.data) .*  ( H*ones(1,size(glm.data,1)) ) ) * ( (plm*ones(1,size(siderealTime,2))) .*timeRotation(glm.Lmax,siderealTime,true) );
CP = transpose( (transpose(timeRotation(glm.Lmax,siderealTime,true)) .* (ones(size(siderealTime,2),1) * transpose(plm)) ) * ( (ones(size(glm.data,1),1) * transpose(H)) .* glm.data ));

if numFreqs~=glm.numFreqs || flow~=glm.flow || deltaF~=glm.deltaF
  % need to interpolate
  freqCP = flow + (0:numFreqs-1)'*deltaF;
  CP=interp1(freq,CP,freqCP);
end

return
