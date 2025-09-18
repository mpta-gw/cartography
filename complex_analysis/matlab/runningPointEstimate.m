function runningPointEstimate(ccStatsFileName, badGPSTimes, ...
                              resampleRate, segmentDuration, ...
                              outputFileNamePrefix, figureLegend, doOverlap, displayResults)
%
%  runningPointEstimate --- generates running point estimate
%
%  runningPointEstimate(ccStatsFileName, badGPSTimes, resampleRate,
%  segmentDuration, outputFileNamePrefix, figureLegend, doOverlap) 
%  generates a plot of the running point estimate over the duration 
%  of the run.
%
%  Input:
%
%    ccStatsFileName = file containing the CC stat values and 
%      theoretical sigmas for the run
%    badGPSTimes = a column vector of 'bad' GPS times corresponding 
%      to outlier segments (to be ignored from the analysis)
%    resampleRate = sample rate (Hz)  (typically, 1024)
%    segmentDuration = segment duration (sec) (typically, 60)
%    outputFileNamePrefix = prefix of filename (and figure) that will 
%      contain the running point estimate values
%    figureLegend = text for figure legend
%    doOverlap = 0 standard weighting by 1/sigma_I^2
%              = 1 combine data using 50% overlapping windows  
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%   
%  $Id: runningPointEstimate.m,v 1.31 2006-09-11 23:14:55 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyyymmddhhMMss = datestr(now,31);
yyyymmdd  = datestr(now,29);

% read in cc stats from file
[gpsTimes, stats, sigmas] = ...
    textread(ccStatsFileName, '%f%f%f\n', -1, 'commentstyle', 'matlab');
if isempty(gpsTimes)
  return;
end;
 
% sort the arrays on gps start times (if not already sorted)
[gpsTimes, ind] = sort(gpsTimes);
stats = stats(ind);
sigmas = sigmas(ind);
                                                                                
% ignore bad gps times
[gpsTimes, ind] = setdiff(gpsTimes, badGPSTimes);
stats = stats(ind);
sigmas = sigmas(ind);

% normalise to get pt estimates and error bars for each segment
pointEsts = stats/segmentDuration;
errorBars = sigmas/segmentDuration;

% check that there's still exists data to combine
if isempty(gpsTimes) return; end
                                                                                
% construct overlapping windows (if necessary)
if doOverlap==1
  window1=hann(resampleRate*segmentDuration);
  window2=hann(resampleRate*segmentDuration);
  [w2bar, w4bar, woverlap4bar] = windowFactors(window1, window2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize cumulative arrays
pointEsts_cum = zeros(length(gpsTimes),1);
errorBars_cum = zeros(length(gpsTimes),1);
pointEsts_cum(1) = pointEsts(1);
errorBars_cum(1) = errorBars(1);

% combine the data for all other segments
ii = 2;
while ii<=length(gpsTimes)

  %fprintf('analysing segment %d\n', ii);

  % initialisation for possible overlapping segments
  jj=1;
  numerator_o = pointEsts(ii-1)/(errorBars(ii-1)^2);
  denominator_o = 1/(errorBars(ii-1)^2);
  pointEst_o = numerator_o/denominator_o;
  errorBar_o = sqrt(1/denominator_o);
  numerator_e = 0;
  denominator_e = 0;
  denominator_oe = 0;                                                                             
  while ( (doOverlap==1) & (gpsTimes(ii)==gpsTimes(ii-1)+segmentDuration/2) ) 
    % save index of last non-overlapping data segment on first
    % pass through (NOTE: might be 0)
    if jj==1
      nonovl_index = ii-2;
    end;

    % increment jj
    jj=jj+1;

    % update info for even and odd overlapping segments
    if mod(jj,2)==0
      numerator_e = numerator_e + ...
                    pointEsts(ii)/(errorBars(ii)^2);
      denominator_e = denominator_e + 1/(errorBars(ii)^2);
      pointEst_e = numerator_e/denominator_e;
      errorBar_e = sqrt(1/denominator_e);
    else
      numerator_o = numerator_o + ...
                    pointEsts(ii)/(errorBars(ii)^2);
      denominator_o = denominator_o + 1/(errorBars(ii)^2);
      pointEst_o = numerator_o/denominator_o;
      errorBar_o = sqrt(1/denominator_o);
    end
    denominator_oe = denominator_oe + ...
                     1/(errorBars(ii-1)^2) + 1/(errorBars(ii)^2);
 
    % calculate components of variance covariance matrix
    C_oo = 1/denominator_o;
    C_ee = 1/denominator_e;
    C_oe = C_oo*C_ee*0.5*(woverlap4bar/w4bar)*0.5*denominator_oe;
    detC = C_oo*C_ee - C_oe*C_oe;
                                                                                
    % construct optimal weighting factors from the covariance matrix
    lambda_o = (C_ee - C_oe)/detC;
    lambda_e = (C_oo - C_oe)/detC;

    % optimal combination for overlapping data
    pointEst_ovl = (pointEst_o*lambda_o + pointEst_e*lambda_e)/... 
                   (lambda_o + lambda_e) ;
    errorBar_ovl = sqrt(1/(lambda_o+lambda_e));                                                             
    % combine with previous non-overlapping cumulative pt estimate
    % as simple weighted sum (provided the index of the previous 
    % non-overlapping point estimate is not zero)
    if nonovl_index==0
      pointEsts_cum(ii) = pointEst_ovl;
      errorBars_cum(ii) = errorBar_ovl;
    else
      numerator = pointEsts_cum(nonovl_index)/ ...
                  (errorBars_cum(nonovl_index)^2) + ...
                  pointEst_ovl/(errorBar_ovl^2);
      denominator = 1/(errorBars_cum(nonovl_index)^2) + 1/(errorBar_ovl^2);
      pointEsts_cum(ii) = numerator/denominator;
      errorBars_cum(ii) = sqrt(1/denominator);
    end

    % increment ii
    ii=ii+1;

    if ii>length(gpsTimes)
      break;
    end

  end % while loop

  % need this to keep the while loop from executing when ii>lenght(times)
  if ii>length(gpsTimes)
    break;
  end

  % non-overlapping data (simple weighted average)
  numerator = pointEsts_cum(ii-1)/(errorBars_cum(ii-1)^2) + ...
              pointEsts(ii)/(errorBars(ii)^2);
  denominator = 1/(errorBars_cum(ii-1)^2) + 1/(errorBars(ii)^2);
  pointEsts_cum(ii) = numerator/denominator;
  errorBars_cum(ii) = sqrt(1/denominator);

  % increment ii
  ii=ii+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display final values
fprintf('Final point estimate  = %e\n', pointEsts_cum(end));
fprintf('Final error bar   = %e\n', errorBars_cum(end));
fprintf('Number segs = %d\n', length(gpsTimes));

% write results to file
filename = [outputFileNamePrefix '_runningPointEstimate.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', yyyymmddhhMMss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\t%s\n', '%gps time', 'cumul pt est', 'error bar');
for ii=1:length(gpsTimes)
  fprintf(fid, '%d\t%e\t%e\n', gpsTimes(ii), pointEsts_cum(ii), errorBars_cum(ii));
end 
fclose(fid);

% Calculate days since start of run
days = (gpsTimes - gpsTimes(1))/86400;

% plot
h = figure;
plot(days, pointEsts_cum, 'k.', ...
     days, pointEsts_cum + 1.65*errorBars_cum, 'g.', ...
     days, pointEsts_cum - 1.65*errorBars_cum, 'b.');
grid on;

% Use this instead of the ylim() fn because ylim() doesn't work in compiled matlab
g = get(h, 'CurrentAxes');
set(g, 'ylim', ...
  [pointEsts_cum(end)-20*errorBars_cum(end) pointEsts_cum(end)+20*errorBars_cum(end)]);
xlabel('Days since start of run');
ylabel('Point estimate +/- 1.65\sigma');
legend(figureLegend);
title(['runningPointEstimate ' yyyymmdd]);
filename = [outputFileNamePrefix '_runningPointEstimate'];
format compact;
print(h, '-depsc2', filename);
print(h, '-dpng', filename);
format;
saveas(h, filename, 'fig');

if (~displayResults)
  close(h);
end;
% Check whether we're running compiled matlab
try
  deployedFlag = isdeployed;
catch
  deployedFlag = 0;
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 2 - runningPointEstIntegrand
integrand = pointEsts./(errorBars.^2);

h = figure;
plot(days,integrand,'b.');
xlabel('Days since start of run');
ylabel('\Omega/\sigma^2');

title(['runningPtEstIntegrand ' yyyymmdd]);
filename = [outputFileNamePrefix '_runningPtEstIntegrand'];
format compact;
print(h, '-depsc2', filename);
print(h, '-dpng', filename);
format;
saveas(h, filename, 'fig');

if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One more plot -- just the errorBars as a function of time.

% Change to units of Days
dataanalyzed = segmentDuration*[1:length(errorBars_cum)]/86400;

h = figure;
semilogy(dataanalyzed,errorBars_cum,'b-');
xlabel('Data analyzed (days)')
ylabel('\sigma')

hold on;
semilogy(dataanalyzed,7.5e-7*sqrt(365./dataanalyzed),'r-');

legend('\sigma','\sigma_{SRD}');
title(['runningSigma ' yyyymmdd]);
filename = [outputFileNamePrefix '_runningSigma'];
format compact;
print(h, '-depsc2', filename);
print(h, '-dpng', filename);
format;
saveas(h, filename, 'fig');

if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

return;
