function [P, X, Fisher, Sky, plm, P1, P2, h1h2] = plmToPlm (sourceFile, lmax)

% arguments:
% 
% lmax           Highest spherical harmonic to use in the search
%
% sourceFile     Real 3 column matrix or one column vector.
%                Power for L=0,...,Lmax
%                and M=-L,..L for f = flow+deltaF*[0:numFreqs-1]
%                packed as a 2-d matrix as shown below:
%
%              |   l   m   power
%   ---------------------------------
%   M=-L,L=Lmax|
%   .......    |
%   M=-1,L=Lmax|
%   .......    |
%   M=-1,L=1   |
%   M=0,L=0    |
%   M=0,L=1    |
%   .......    |
%   M=0,L=Lmax |
%   M=1,L=1    |
%   .......    |
%   M=1,L=Lmax |
%   .......    |
%   M=L=Lmax   
%
%           or, just the 3rd column of that matrix
%
%      *** make sure your source is a real spherical harmonic ***
%                real P(Omega)P*(l,m) = (-1)^m P(l,-m)
%
%  output:   
%	
%       P         - Power of each spherical harmonic, in same form as
%                   calGammaLM.m
%       X         - Output of calX.m
%       Fisher    - Output of calFisher.m
%       Sky       - a struct containing the output of doSpH.m
%       plm       - Power of each spherical harmonic, in real spherical
%                   harmonics (lmcosi form)
%
%  Routine written by Madeleine Udell.
%  Contact madeleine_udell@yale.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% First we set up a lot of parameters %%%%%%%%%%


gammaLM_coeffsPath='/archive/home/sballmer/stochastic/src/matapps/src/searches/stochastic/CrossCorr';
% Time
%Assuming integration for one day and 192sec chunk
% Assuming one sidereal day is 86400s, instead of 86164s
totalTime = 86164.0;
deltaTime = 192.0*1;
timeBins = deltaTime/2:deltaTime:totalTime;
timeBins = (0:1:23)*3600;

% Frequency
% Assuming 512Hz upper cut frequency and 1/4Hz bin size
upperFreq=512.0;
deltaFreq=0.25;
freqBins = deltaFreq/2:deltaFreq:upperFreq;
flow=freqBins(1);
numFreqs=length(freqBins);
deltaF=deltaFreq;

% Frequency spectrum of sources (power law)
alpha=0;
HAlpha=1;
detectorPair='HL';
if strcmp(detectorPair,'HL');
    det1='LHO';
    det2='LLO';
    dubdet1='HH';
    dubdet2='LL';
else error('Unrecognized Detector Pair');
end;

%%% parameters for testSpH.m %%%
w1w2bar=1;
w1w2squaredbar=1;
FreqIntFlag=true;
outputFilePrefix='/archive/home/sballmer/stochastic/sgwb/S5/output/S5_SH_1job_Y2Y';
maxSegmentsPerFile=60;
jobNr=1;
trialNr=1;
try, lmax; catch lmax=3; end;

lsource=sqrt(length(sourceFile))-1;

% The Power Spectra are currently set to one for the whole frequency range
% and at all times.
% You might want to make them something more realistic


%%%%%%%%%% Now we start in on the chain of programs called %%%%%%%%%%

% injectSpH.m
% outputs 
% H          Real vector. Frequency power spectrum of the source
% invH       Real vector. Inverse of frequency power spectrum of the source
% h1h2       Complex matrix. Contains h1*(f)h2(f) at every time bin
%            Index 1: time, Index 2: frequency
% table      Should be the same as map

% We add deltaTime/2 because the doSpH analysis evaluates at the middle of the time intervals
siderealTime=GPStoGreenwichMeanSiderealTime(timeBins+deltaTime/2);
  
 freq = flow + (0:numFreqs-1)'*deltaF;
 H = HAlpha*((freq/100.0) .^ alpha);
 glm = calGammaLM(gammaLM_coeffsPath,detectorPair, lsource, numFreqs, flow, deltaF);
 glm1= calGammaLM(gammaLM_coeffsPath,dubdet1,      lsource, numFreqs, flow, deltaF);
 glm2= calGammaLM(gammaLM_coeffsPath,dubdet2,      lsource, numFreqs, flow, deltaF); 

 h1h2   = SpH2CrossPower(glm, sourceFile, siderealTime, numFreqs, flow, deltaF, HAlpha, alpha);
 P1     = SpH2CrossPower(glm1,sourceFile, siderealTime, numFreqs, flow, deltaF, HAlpha, alpha);
 P2     = SpH2CrossPower(glm2,sourceFile, siderealTime, numFreqs, flow, deltaF, HAlpha, alpha);

P1=ones(size(P1));%+randn(size(P1)).^2;
P2=ones(size(P1));%+randn(size(P2)).^2;


% vSpH.m
% outputs a struct vSpH
tic; fprintf('\nConstructing vSpH ');
vSpH=SpH(gammaLM_coeffsPath,detectorPair,lmax,numFreqs,flow,deltaFreq,H,w1w2bar,w1w2squaredbar,deltaTime,FreqIntFlag,outputFilePrefix,maxSegmentsPerFile,jobNr,trialNr);
fprintf('DONE\n'); toc;

% doSpH.m
% outputs a struct
tic;
for ii=1:length(timeBins)
  GPSstart=timeBins(ii);
  fprintf('\nCalling doSpH for t=%d\n',GPSstart);
  vSpH=doSpH(vSpH,h1h2(:,ii),P1(:,ii),P2(:,ii),GPSstart);
end
toc;

% saveSpHSet.m
% outputs a struct vSpH and a file
vSpH=saveSpHSet(vSpH);

% invertX.m
% sums the Fisher matrices and the X vectors from each time segment and
% uses them to calculate the coefficients of the spherical harmonics.
% saves them as a struct in a file outputFilePrefix_SpH.job1.trial1.mat
indexfilename=sprintf('%s_SpH.job1.trial1.mat',outputFilePrefix);
tic; fprintf('\nInverting Fisher matrix and saving output in %s ', indexfilename);
Sky=load(indexfilename);
[P, Fisher, X, invFisher] = invertX(indexfilename);
%save(indexfilename,'Sky', 'X', 'Fisher', 'P');
fprintf('DONE\n');toc;

plm=plm2plmreal(P);
