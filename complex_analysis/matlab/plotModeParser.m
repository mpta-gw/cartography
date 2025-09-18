function mode = plotModeParser(modestr)
% 
% this function is used to translate strings into 
% a plot mode for plotMapAitoff and plotPlmMap
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

if isstr(modestr)
    modestr=upper(modestr);
    cl =strfind(modestr,'CL');
    cal=strfind(modestr,'CAL');
    if strfind(modestr,'ALL')
        mode=0;
    elseif strfind(modestr,'SNR')
        mode=1;
    elseif strfind(modestr,'PE')
        mode=2;
    elseif strfind(modestr,'Y')
        mode=2;
    elseif strfind(modestr,'SIG')
        mode=3;
    elseif strfind(modestr,'LOG')
        mode=4;
    elseif cl
        if cal
            mmm=str2num(modestr(cal+3:end));
            nnn=str2num(modestr(cl+2:cal-1));
            if nnn<=0.5 || nnn>=1.0
                error('Confidence level has to be between 0.5 and 1');
            end
            if mmm<=0.0 || mmm>=1.0
                error('Calibration uncertainty has to be between 0 and 1');
            end
            mode=round(nnn*1000)+mmm;
        else
            nnn=str2num(modestr(cl+2:end));
            if nnn<=0.5 || nnn>=1.0
                error('Confidence level has to be between 0.5 and 1');
            end
            mode=round(nnn*100);
        end
    else
        mode=-abs(str2num(modestr));
    end
else
    mode=modestr;
end