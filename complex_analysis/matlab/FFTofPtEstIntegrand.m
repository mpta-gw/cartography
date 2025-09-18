function [t, omega_t] = FFTofPtEstIntegrand(ptEstIntegrand, ...
                                            resampleRate, ...
                                            outputFileNamePrefix, ...
                                            displayResults, ...
                                            figureNumber, figureLegend, ...
                                            complexflag);
%
%  [t, omega_t] = FFTofPtEstIntegrand(ptEstIntegrand, resampleRate, ...
%                                     outputFileNamePrefix, ...
%                                     displayResults, ...
%                                     figureNumber, figureLegend, ...
%                                     complexflag)
%
%  calculates and plots (if desired) the FFT of the integrand of the
%  optimal point estimate.  
%
%  Input:
%
%    ptEstIntegrand = a freq-series data structure containing the
%      integrand of the point estimate of Omega
%    resampleRate = resample rate of the data used to generate the
%      integrand of the point estimate
%    outputFileNamePrefix = prefix of filename that will contain
%      the FFT of the integrand of the point estimate; 
%      prefix of filenames for figures if displaying results.
%    displayResults = 0 don't display results to screen (default = 1)
%    figureNumber = number for figure (if displayResults = 1)
%    figureLegend = text for figure legend (if displayResults = 1)
%    complexFlag = if true, sum only over positive frequencies
%                  if false, infer negative frequencies to ensure reality
%
%  Output:
%
%    t = discrete times (-1/(2*deltaF) ... 1/(2*deltaF) in steps of
%      1/resampleRate  (NOTE: deltaF is the frequency resolution of 
%      the integrand of the point estimate) 
%    omega_t = time-series containing the FFT of integrand of the point 
%      estimate of Omega
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk or john.whelan@ligo.org
%  
%  $Id: FFTofPtEstIntegrand.m,v 1.21 2006-08-18 01:51:40 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check complex flag
try
  complexflag;
catch
  complexflag = false;
end;

yyyymmddhhMMss  = datestr(now,31);
yyyymmdd = datestr(now,29);

% extract freq series data
[data, flow, deltaF] = extractFreqSeries(ptEstIntegrand);

numFreqs = length(data);
f = flow+deltaF*transpose([0:numFreqs-1]);
fhigh = flow+deltaF*numFreqs;

% some derived quantities
deltaT = 1/resampleRate;
fNyq = resampleRate/2;

% extend frequency range from -fNyq to +fNyq
numFreqs_pre = floor(flow/deltaF)-1;
f_pre = deltaF*[1:numFreqs_pre]';
numFreqs_post = floor((fNyq - fhigh)/deltaF);
f_post = fhigh + deltaF*[0:numFreqs_post-1]';
fp =  [f_pre; f; f_post];
fn = -flipud(fp);
f_tot = [fn; 0; fp];

% zero pad cc spectrum, extending to negative frequencies
integrand_pre  = zeros(numFreqs_pre,1);
integrand_post = zeros(numFreqs_post,1);
integrand_p = [integrand_pre; data; integrand_post];

% infer values at neg freqs unless statistics complex
if (complexflag == false)
  integrand_n = flipud(conj(integrand_p));
else
  integrand_n = zeros(size(integrand_p));
end;

integrand_tot = [0; integrand_p; integrand_n];

% take FFT
fft_integrand = fftshift(fft(deltaF*integrand_tot)); % NOTE: deltaF
%fprintf('maximum real part = %e\n',max(abs(real(fft_integrand))) );
%fprintf('maximum imag part = %e\n',max(abs(imag(fft_integrand))) );
if (complexflag == false)
  fft_integrand = real(fft_integrand);  % NOTE: fft_integrand is real to
                                        % round-off if complexflag
                                        % is false, so i'm taking the
                                        % real part ONLY so that matlab
                                        % doesn't complain later on 
                                        % when making plots
end;
  
omega_t=flipud(fft_integrand); % NOTE: flipud for timeshift convention

t = -1/(2*deltaF) + deltaT : deltaT : 1/(2*deltaF) -deltaT;
t = transpose(t); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write FFT of point estimate integrand to file
                                                                                
filename = [outputFileNamePrefix '_FFTofPtEstIntegrand.dat'];
fid = fopen(filename, 'w');
fprintf(fid, '%% Date and time of this run: %s\n', yyyymmddhhMMss);
fprintf(fid, '%s\n', '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%s\t%s\n', '%time lag (s)', 'FFT of pt est integrand');
for ii=1:length(t)
  fprintf(fid, '%d\t%e\n', t(ii), omega_t(ii));
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display results (if desired)

% Don't try this with complex statistics for now
if complexflag
  return
end;

% value at zero-lag
ind = find(t==0);
pointEstimate = omega_t(ind);
fprintf('Pt estimate (zero lag) = %e\n', pointEstimate);

% plot
h = figure(figureNumber);
plot(t, omega_t);
grid on;
xlabel('Lag (sec)','FontSize',14)
ylabel('FFT of Integrand of Pt Estimate','FontSize',14)
title(['FFTofPtEstIntegrand ' yyyymmdd]);
legend(figureLegend);
filename = [outputFileNamePrefix '_FFTofPtEstIntegrand'];
format compact;
print(h, '-depsc2', filename);
print(h, '-dpng', filename);
saveas(h, filename, 'fig');
format;
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

return
