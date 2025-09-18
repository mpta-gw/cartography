function ...
  [sky, numSegments] = ...
    combineResultsRM(jobCell, ...
                      segmentDuration, first, last, ...
                      doOverlap, window1, window2)
%
%  
%  function ...
%    [sky, numSegments] = ...
%      combineResultsRM(jobCell, ...
%                        segmentDuration, first, last, ...
%                        doOverlap, window1, window2)
%
%  optimally combines cc stats, theor sigmas, cc spectra, and 
%  sensitivity integrands.
%
%  Input:
%
%    jobCell - cell array (1 entry per segment)
%              each entry is a struct
%               - time: GPS time
%               - data: Nx2 array, 1st column = ccStat, 2nd column = sigma
%                       rows correspond to different points in the sky
%              instead of data it can also contain filename and segmentOffset:
%               - filename      = 'file with the actual data'
%               - segmentOffset = index offset
%    segmentDuration = length of analysis segment in sec (typically 60)
%    first, last - first and last segment inder in jobCell
%
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows
%    window1 = array containing the window used for the first time-series
%              (should have numPoints = segmentDuration*resampleRate)
%    window2 = array containing the window for the second time-series
%              (should have numPoints = segmentDuration*resampleRate)
%
%  Output:
%
%    sky      = Nx2 array with optimal point estimate of Omega0 (column 1) and 
%              theoretical error bar for the point estimate (column 2)
%              rows correspond to different points in the sky
%    numSegments = number of data segments combined
%
%  Routine written by Joseph D. Romano and Stefan Ballmer
%  Contact Joseph.Romano@astro.cf.ac.uk / sballmer@ligo.mit.edu
%
%  $Id: combineResultsRM.m,v 1.3 2006-05-08 18:02:15 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract number of segments
numSegments = last-first+1;
loadedFile='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  seg=first;
  % get the appropriate data
  try data=jobCell{seg}.data; catch
    if ~strcmp(loadedFile,jobCell{seg}.filename)
      loadedFile=jobCell{seg}.filename; clear p;
      try p=load(loadedFile); catch fprintf('Could not load file %s\n',loadedFile); sky=[]; numSegments=0; return; end;
    end
    data=p.Sky{seg-jobCell{seg}.segmentOffset}.data;
  end
  N=size(data,1);

% apply optimal weighting
if ( (doOverlap == false) | (numSegments == 1) )
  Y=zeros(N,1);
  E=zeros(N,1);
  for seg=first:last
    % get the appropriate data
    try data=jobCell{seg}.data; catch
      if ~strcmp(loadedFile,jobCell{seg}.filename)
        loadedFile=jobCell{seg}.filename; clear p;
        try p=load(loadedFile); catch fprintf('Could not load file %s\n',loadedFile); sky=[]; numSegments=0; return; end;
      end
      data=p.Sky{seg-jobCell{seg}.segmentOffset}.data;
    end
    if length(data)>0
      invVar=1./(data(:,2).^2);
      if N==0
        N=size(data,1);
	Y=zeros(N,1);
        E=zeros(N,1);
      end
      Y=Y+data(:,1).*invVar;
      E=E+invVar;
    end;
  end
  sky=[Y./E,1./sqrt(E)]/segmentDuration;
  clear invVar E Y;
else
  % even and odd point estimates
  Y_o=zeros(N,1);
  Y_e=zeros(N,1);
  E_o=zeros(N,1);
  E_e=zeros(N,1);
  sigmaratio=zeros(N,1);
  [w2bar, w4bar, woverlap4bar] = windowFactors(window1, window2);
  wfac=0.5*(woverlap4bar/w4bar)*0.5;
  for seg=first:last
    % get the appropriate data
    try data=jobCell{seg}.data; catch
      if ~strcmp(loadedFile,jobCell{seg}.filename)
        loadedFile=jobCell{seg}.filename; clear p;
        try p=load(loadedFile); catch fprintf('Could not load file %s\n',loadedFile); sky=[]; numSegments=0; return; end;
      end
      data=p.Sky{seg-jobCell{seg}.segmentOffset}.data;
    end
    invVar=1./(data(:,2).^2);
    if mod(seg-first,2)==0
      Y_o=Y_o+data(:,1).*invVar;
      E_o=E_o+invVar;
    else
      Y_e=Y_e+data(:,1).*invVar;
      E_e=E_e+invVar;
    end
    sigma2J  = data(:,2).^2;
    if seg ~= first
      sigma2IJ = wfac * (sigma2I+sigma2J);
      sigmaratio = sigmaratio + sigma2IJ./(sigma2I.*sigma2J);
    end
    sigma2I=sigma2J;
  end
  ptEstimate_o = Y_o./E_o/segmentDuration;
  ptEstimate_e = Y_e./E_e/segmentDuration;
  % construct covariance matrix for the even and odd data
  C_oo = (1/segmentDuration^2)./E_o;
  C_ee = (1/segmentDuration^2)./E_e;
  sigmaratio = sigmaratio * segmentDuration;
  clear invVar E_o Y_o E_e Y_e;  
  C_oe = C_oo.*C_ee.*sigmaratio;
  detC = C_oo.*C_ee - C_oe.*C_oe;

  % construct optimal weighting factors from the covariance matrix
  lambda_o = (C_ee - C_oe)./detC;
  lambda_e = (C_oo - C_oe)./detC;
  
  % final point estimate
  % final error bar 
  variance = 1./(lambda_o + lambda_e);
  sky = [(ptEstimate_o.*lambda_o + ptEstimate_e.*lambda_e) .* variance , sqrt(variance) ];
end
return

