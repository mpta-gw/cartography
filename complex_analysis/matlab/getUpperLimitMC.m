function UL = getUpperLimitMC(data, Y, sigma, confidence, calError)
% function UL = getUpperLimitMC(data, Y, sigma, confidence, calError)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fuction is an extension to getUpperLimit.m,
% i.e. it marginalizes over the calibration error
% and calculates the upper limit.
% But instead of assuming a Gaussian distribution,
% it take a sample data (trials) as first argument.
%
% Input arguments:
%  data      : array [Ntrial x n] or [Ntrial x 1]
%              n is the number of [Y,sigma] pairs
%  Y         : point estimate [n x 1]
%  sigma     : error bar [n x 1]
%  confidence: confidence level (e.g. 0.9 for 90%)
%  calError  : calibration error (e.g. 0.1 for 10%)
%              currently limited to <0.2 (due to integrand singularity)
%
% Output argument
%  UL        : upper limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eric Thrane & Stefan Ballmer 2011 Feb 23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gaussFlag=or(~length(data),ischar(data));

% adjust the data size
nn=length(Y);
if and(nn>1,size(data,2)==1)
    data=data*ones(1,nn);
end

L = size(Y,1);
for n=1:L
  % Define the domain over which the likelihood function is defined
  Omg = 0 : sigma(n)/1000 : 10*sigma(n);

  % testing and debugging flag
  Marginalize = true;
%  Marginalize = false;

  if Marginalize
    % marginalized likelihood
    LklhdM = zeros(size(Omg));
    for lambda=1-5*calError:calError/50:1+5*calError
      % likelihood function defined such that the mean is Y
      %Lklhd0 = hist(data(:,n)+Y(n), Omg);
      %Lklhd0=Lklhd0/sum(Lklhd0);

      % likelihood function given lambda
      % ... but now scale mean and sigma by lambda
      if gaussFlag
          Lklhd = exp(-1/2*(Y(n)-lambda*Omg).^2/sigma(n)^2);
      else
          Lklhd = hist( (data(:,n)+Y(n))/lambda, Omg)./lambda;
          % kill the overflow
          Lklhd(  1)=0;
          Lklhd(end)=0;
      end

%WRONG:      % noramlize likelihood function given lambda
%WRONG:      Lklhd = Lklhd/sum(Lklhd);
      

      % define prior-----------------------------------------------------------
      % Gaussian prior (may not formally converge)
      prior=exp(-0.5*(lambda-1)^2/calError^2);

      % flat prior
      % prior = (heaviside(lambda-1+calError).*heaviside(calError-lambda+1)) ...
      %  /(2*calError);
      % heaviside functions produce NaNs at heaviside=0
      % the solution is to make it a gradual step
      % if(isnan(prior)), prior=1/(4*calError); end;
      %------------------------------------------------------------------------

      % calculate marginalized likelihood
      LklhdM = LklhdM + Lklhd * prior;
      if any(isnan(Lklhd))
          disp('NAN');
      end

      % How to calculate mean and sigma of Lklhd distribution.
      %meanL0 = sum(Omg.*Lklhd0);
      %meanL = sum(Omg.*Lklhd);
      %stdL0 = sqrt( sum(Omg.^2.*Lklhd0) - meanL0^2 );
      %stdL = sqrt( sum(Omg.^2.*Lklhd) - meanL^2 );
    end
  else
    % no marginalization
    % likelihood function defined such that the mean is Y
    Lklhd0 = hist(data(:,n)+Y(n), Omg);
    LklhdM = Lklhd0;
  end

  % calculate posterior assuming positivity
  posterior = LklhdM;
  posterior(Omg<0) = 0;

  % define cumulative distribution
  cdf = cumsum(posterior,2);
  % normalize
  cdf=cdf./cdf(end);

  % find index where cdf crosses confidence
  idx = abs(cdf-confidence)==min(abs(cdf-confidence));

  if sum(idx)>1
    UL(n) = sum(Omg(idx))/sum(idx);
    fprintf('test message: %i\n', sum(idx));
  else
    UL(n) = Omg(idx);
    %indx=find(idx);
    %UL(n)=interp1(cdf(indx-1:indx+1),Omg(indx-1:indx+1),confidence);
  end
end

return
