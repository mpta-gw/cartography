function [o_map,o_sigmaPix]=plotPlmMap(plm,covar,mode,scale,RadiometerFlag);
% plots a spherical harmonics map
% The function supports different syntaxes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) only 1 argument: plm
%    --> A map of plm is plotted (backward compatibility mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B) 2 arguments given: plm, cover (i.e. the covariance matrix of plm)
%    --> The SNR map is plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C) 3 arguments given: plm, covar, mode
%    --> plot depends on the value of mode:
%
%        mode   |   action
% ------------------------------------------------------------------------
%         0     |   plot SNR, Point-Estimate and Sigma map (mode 1,2,3)
%         1     |   plot an SNR map (default if mode is not specified)
%         2     |   plot Point-Estimate map
%               |   similar to syntax A), but adds a title
%         3     |   plot Sigma map
%         4     |   plot logLikelihood '-log(p(x>abs(snr)))'
%   50<mode<100 |   plot bayesian UL with 'confidence=mode/100'
%               |   e.g. mode=90 is a 90% C.L. UL map
%  500<mode<1000|   plot bayesian UL with 'confidence=floor(mode/1000)'
%               |   and calibration uncertainty given behind the decimal
%               |   point e.g. mode=900.12 is a 90% C.L. UL map for
%               |   12% calibration uncertainty
%         -1    |   plot Point-Estimate map, same as syntax A)
%         -2    |   plot Sigma map, but no title
%
%    Note 1: The mode argument is passed to plotMapAitoff and parsed there.
%    Note 2: instead of mode numbers, simple strings can be used.
%            see plotModeParser.m or below for details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D) 4 arguments: as C), but first apply a scale factor to both plm and
%    covar. This avoid the plotting float overflow that might occur
%    for extremely large or small numbers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E) 5 arguments: as D), 5th arugment is RadiometerFlag
%    RadiometerFlag=true assumes that:
%       plm   = X      (straight from the analysis)
%       covar = Fisher (straight from the analysis)
%    It then plots the requested maps for the radiometer result.
%    All that changes is replacing the call
%       plotMapAitoff([map,sigmaPix],360,181,mode);
%    with a call
%       plotMapAitoff([map./sigmaPix.^2,1./sigmaPix],360,181,mode);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accepted mode strings:
%
% Str  | mode   |   action
% ------------------------------------------------------------------------
% 'ALL'|  0     |   plot SNR, Point-Estimate and Sigma map (mode 1,2,3)
% 'SNR'|  1     |   plot an SNR map (default if mode is not specified)
% 'PE' |  2     |   plot Point-Estimate map
%      |        |   similar to syntax A), but adds a title
% 'SIG'|  3     |   plot Sigma map
% 'LOG'|  4     |   plot logLikelihood '-log(p(x>abs(snr)))'
% 'n'  |  -n    |   plot the n-th map, no title
%      |--------|
% 'CLnnn'       |   plot bayesian UL with 'confidence=mode/100'
%               |   e.g. mode=90 is a 90% C.L. UL map
% 'CLnnn CALmmm'|   plot bayesian UL with 'confidence=floor(mode/1000)'
%               |   and calibration uncertainty given behind the decimal
%               |   point e.g. mode=900.12 is a 90% C.L. UL map for
%               |   12% calibration uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




s=size(plm);

if s(2)==1 && mod(sqrt(s(1)),1)==0
  fprintf('\nInput is vector of complex spherical harmonics coefficients.\n\n');
elseif (s(2)==2 || s(2)==4 ) && mod((sqrt(8*s(1)+1)-3)/2,1)==0
  fprintf('\nInput is vector of real spherical harmonics coefficients.\n\n');
  plm=plmreal2plm(plm);
else
  error('No idea what the input is...');
end

if exist('scale','var')
    plm=plm*scale;
    covar=covar*(scale^2);
end

if not(exist('RadiometerFlag','var'))
    RadiometerFlag=false;
end

if nargin==1
    [map,ra,decl]=makemap(plm,1);
    if nargout==0
        plotMapAitoff(map,360,181,-1);
    else
        o_map=map;
    end
else
    [sigma,sigmaPix]=getSigmaMap(covar,1);
    [map,ra,decl]=makemap(plm,1);
    if ~exist('mode','var'), mode=1; else mode=plotModeParser(mode); end;
    for mm=mode
        if mm==1 || mm==0, fprintf('Plotting SNR map\n'); end;
        if mm==2 || mm==0, fprintf('Plotting Point-Estimate map\n'); end;
        if mm==3 || mm==0, fprintf('Plotting Sigma map\n'); end;
        if mm==-1, fprintf('Plotting Point-Estimate map (no title)\n');end;
        if mm==-2, fprintf('Plotting Sigma map (no title)\n');end;
    end
    if nargout==0
        if RadiometerFlag
            fprintf('Plotting Radiometer map...\n');
            plotMapAitoff([map./(sigmaPix.^2),1./sigmaPix],360,181,mode);
        else
            plotMapAitoff([map,sigmaPix],360,181,mode);
        end
    else
        o_map=map;
        o_sigmaPix=sigmaPix;
    end
end;
