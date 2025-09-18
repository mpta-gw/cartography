gammaLM_coeffsPath='/archive/home/sballmer/stochastic/src/matapps/src/searches/stochastic/CrossCorr';

detPair='HL';
det1=getdetector(detPair(1));
det2=getdetector(detPair(2));
lmax=30;
numFreqs=501;
flow=10;
deltaF=1;
freq = flow + (0:numFreqs-1)'*deltaF;
beta=0;
HBeta=1;

RA=6;
DECL=45;
siderealTime=15;

deg=5;

if 1
  pNP=deltaAtNorthPole(lmax);
  wRot=WignerRotation(lmax,(90-DECL)/180*pi);
  Z=zRotation(lmax,RA/12*pi);
  s=Z*wRot*pNP;
else
  s=randPlm(lmax);
end

gamma=orfIntegrandSymbolic(det1,det2,siderealTime,RA,DECL);
H = HBeta*((freq/100.0) .^ beta);
gammaPt=gamma.gamma0.*exp(i*2*pi*freq.*gamma.tau).*H;


glm = calGammaLM(gammaLM_coeffsPath,detPair, lmax, numFreqs, flow, deltaF);

  [map,RAs,DECLs,dOmega] = makemap(s,deg);
  RADECL=[RAs,DECLs];

CP   = SpH2CrossPower(glm, s,                          siderealTime, numFreqs, flow, deltaF, HBeta, beta);
mCP  = Map2CrossPower(det1, det2, map, dOmega, RADECL, siderealTime, numFreqs, flow, deltaF, HBeta, beta);


close all;
plotPlmMap(s)
figure;
plot(freq,real(gammaPt),'r'  ,freq,real(CP),'b-.',freq,real(mCP),'c:',...
     freq,imag(gammaPt),'m'  ,freq,imag(CP),'k-.',freq,imag(mCP),'g:');
grid
xlabel('Hz');
legend('re gammaPt','re CP','re mCP','im gammaPt','im CP','im mCP',3);

%norm(transpose(h1h2)-CP)./norm(CP)
norm(mCP-CP)./norm(CP)


% reduce lmax
lmaxlow=15;
deg=5;
s2=randPlm(lmaxlow);
disp('starting PlmToPlm');
[P2, X2, Fisher2, Sky2, plm2, P12, P22, h1h22] = plmToPlm    (s2,lmaxlow);
disp('starting PlmToPlmPix');
[P3, X3, Fisher3, Sky3, plm3, P13, P23, h1h23] = plmToPlmPix (s2,lmaxlow,deg);
norm(P2-s2)./norm(s2)
norm(P3-s2)./norm(s2)
norm(P2-P3)./norm(s2)


norm(mCP-CP)./norm(CP)

