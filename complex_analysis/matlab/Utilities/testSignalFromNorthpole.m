
gammaLM_coeffsPath='';

detPair='HL';
lmax=30;
numFreqs=501;
flow=10;
deltaF=1;

glm = calGammaLM(gammaLM_coeffsPath,detPair, lmax, numFreqs, flow, deltaF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal at the north pole

% for P(Omega) = delta fuction at north pole:
pNP=deltaAtNorthPole(lmax);

gammaNP1=transpose(pNP)*glm.data;
freq=flow + (0:numFreqs-1)*deltaF;

det1=getdetector('LHO');
det2=getdetector('LLO');
gamma=orfIntegrandSymbolic(det1,det2,0,0,90);

gammaNP2=gamma.gamma0.*exp(i*2*pi*freq.*gamma.tau);


plotPlmMap(pNP);
figure;
%plot(freq,imag(gammaNP1),freq,ones(size(freq)).*imag(gammaNP2));
plot(freq,abs(gammaNP1-gammaNP2)./abs(gamma.gamma0));
xlabel('Frequency [Hz]');
ylabel('abs');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal at RA=0, DECL=0

R=WignerRotation(lmax,pi/2);

pZP=R*pNP;
gammaZP1=transpose(pZP)*glm.data;

gamma=orfIntegrandSymbolic(det1,det2,0,0,0);
gammaZP2=gamma.gamma0.*exp(i*2*pi*freq.*gamma.tau);

plotPlmMap(pZP);
figure;
%plot(freq,imag(gammaZP1),freq,ones(size(freq)).*imag(gammaZP2));
plot(freq,abs(gammaZP1-gammaZP2)./abs(gamma.gamma0));
xlabel('Frequency [Hz]');
ylabel('abs');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And a random direction
RA=135 /15; % in h
DECL=-22.5; % in deg

RR=WignerRotation(lmax,(90-DECL)/180*pi);

ZZ=zRotation(lmax,RA/12*pi);
pRN=ZZ*RR*pNP;
gammaRN1=transpose(pRN)*glm.data;

gamma=orfIntegrandSymbolic(det1,det2,0,RA,DECL);
gammaRN2=gamma.gamma0.*exp(i*2*pi*freq.*gamma.tau);

plotPlmMap(pRN);
figure;
%plot(freq,imag(gammaRN1),freq,ones(size(freq)).*imag(gammaRN2));
plot(freq,abs(gammaRN1-gammaRN2)./abs(gamma.gamma0));
xlabel('Frequency [Hz]');
ylabel('abs');
grid on


