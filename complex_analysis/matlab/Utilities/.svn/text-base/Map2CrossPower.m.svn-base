function CP = Map2CrossPower(det1, det2, map, dOmega, RADECL, siderealTime, numFreqs, flow, deltaF, H, beta)

% calculates the cross-power due to a background given by H*map
% input: det1, det2   detector1 and detector 2 (use getdetector.m)
%        map          vector of P(Omega)
%        dOmega       vector of angular weights (can be scalar, e.g 1)
%        RADECL       right ascension in hours and declination in degrees [Nx2 vector]
%        siderealTime evaluate at this sidereal time
%        numFreqs     number of frequency bins
%        flow         lowest frequency
%        deltaF       frequency bin separation
%        H            power spectrum, either vector with numFreq elements,
%                     or value at 100Hz (scalar)
%        beta         power law (if H is scalar), default =0
%
% output:CP           Cross-Power

% optimized for speed

try, beta; catch beta=0; end
numTimes=length(siderealTime);

if numTimes==1
  freq = flow + (0:numFreqs-1)'*deltaF;
  if length(H)==1
    H=H*((freq/100.0) .^ beta);
  end

  mapdO=map.*dOmega;
  M=size(RADECL,1);
  CP=zeros(numFreqs,1);
  s = det2.r - det1.r;
  % distance between sites
  distance = norm(s);
  % unit vector Omega
  c=299792458;
  w=pi/12;
  for ii=1:M
    psi=w*(siderealTime-RADECL(ii,1));
    theta=-pi/2+pi/180*RADECL(ii,2);
    ctheta=cos(theta);                stheta=sin(theta);
    cpsi  =cos(psi);                  spsi  =sin(psi);
    Omega1=-cpsi*stheta;              Omega2=spsi*stheta;                 Omega3=ctheta;
    pp11=ctheta^2*cpsi^2-spsi^2;      pp12=-(ctheta^2+1)*cpsi*spsi;       pp13=ctheta*cpsi*stheta;
                                      pp22=ctheta^2*spsi^2-cpsi^2;        pp23=-ctheta*spsi*stheta;
                                                                          pp33=stheta^2;
    px11=-2*ctheta*cpsi*spsi;         px12=ctheta*(spsi^2-cpsi^2);        px13=-stheta*spsi;
                                      px22=2*ctheta*cpsi*spsi;            px23=-stheta*cpsi;
                                                                         %px33=0;
    g.F1p=det1.d(1,1)*pp11+2*det1.d(1,2)*pp12+2*det1.d(1,3)*pp13+det1.d(2,2)*pp22+2*det1.d(2,3)*pp23+det1.d(3,3)*pp33;
    g.F2p=det2.d(1,1)*pp11+2*det2.d(1,2)*pp12+2*det2.d(1,3)*pp13+det2.d(2,2)*pp22+2*det2.d(2,3)*pp23+det2.d(3,3)*pp33;
    g.F1x=det1.d(1,1)*px11+2*det1.d(1,2)*px12+2*det1.d(1,3)*px13+det1.d(2,2)*px22+2*det1.d(2,3)*px23;%+det1.d(3,3)*px33;
    g.F2x=det2.d(1,1)*px11+2*det2.d(1,2)*px12+2*det2.d(1,3)*px13+det2.d(2,2)*px22+2*det2.d(2,3)*px23;%+det2.d(3,3)*px33;
    % divide by 2 since Hf is defined as sum over pol of power
    g.gamma0=(g.F1p * g.F2p + g.F1x * g.F2x)/2;	     
    g.tau=(Omega1*s(1)+Omega2*s(2)+Omega3*s(3))/c;
  
    g.data=g.gamma0*exp(i*2*pi*freq*g.tau);
    CP=CP+mapdO(ii)*g.data;
  
    % lines for testing against slower code
    %gt=orfIntegrandSymbolic(det1,det2,siderealTime,RADECL(ii,1),RADECL(ii,2));
    %v=[norm(g.tau-gt.tau),norm(g.F1p-gt.F1p),norm(g.F1x-gt.F1x),norm(g.F2p-gt.F2p),norm(g.F2x-gt.F2x)];
    %fprintf('%d\t%d\t%d\t%d\t%d\n',v);
  end
  CP=CP.*H;


else
  CP=zeros(numFreqs,numTimes);
  for kk=1:numTimes
    CP(:,kk)=Map2CrossPower(det1, det2, map, dOmega, RADECL, siderealTime(kk), numFreqs, flow, deltaF, H, beta);
  end
end
