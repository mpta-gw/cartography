function varargout=iwt(a,d,an,dn,tipe,nvm,pph,intel)
% [x,xr,ts]=iwt(a,d,an,dn,tipe,nvm,pph,intel)
% 
% Performs a wavelet reconstruction.
% based on coefficients 'a' and 'd{n}'.
%
% INPUT :
%
% a      approximation coefficients (after n lowpass)
% d      details (after each of n highpasses)
% tipe   'Daubechies' or 'CDF'
% nvm    Number of vanishing (primal & dual) moments 
%
% OUTPUT:
%
% x      Cell array with reconstruction at different scales;
%        The inverse wavelet transform is given by:
% xr     sum([x{:}],2)
% ts     Time it took to do this.
%
% Examples:
%
% iwt('demo1') through iwt('demo6')
%
% Last modified by fjsimons-at-alum.mit.edu, 11/12/2002

if ~isstr(a)
  t0=clock;
  
  % Default is the CDF 2/4 construction
  defval('tipe','CDF')
  defval('nvm',[2 4])
  defval('pph',3)
  defval('intel',0)
  [h0,f0,P,U,Kp,Ku]=wc(tipe,nvm);
  
  n=length(d);
  
  % Reconstruct the approximation coefficients
  x{1}=sbank(a,h0,f0,P,U,Kp,Ku,'a',n,an,dn,pph,intel);
  % Reconstruct the detail coefficients
  for index=1:n
    x{n+1-index+1}=sbank(d{index},h0,f0,P,U,Kp,Ku,'d',index,...
			 an,dn,pph,intel);
  end
  xr=sum([x{:}],2);
  ts=etime(clock,t0);
  disp(sprintf('Synthesis took %8.4f s',ts))
  varnames={'x' 'xr' 'ts'};
  for index=1:nargout
    varargout{index}=eval(varnames{index});
  end
else
  load('noisdopp'); x0=noisdopp(:);
  switch a
   case 'demo1'
    x0=real(dopnoise(500,200,60,10,70,128));
    [a,d,an,dn]=wt(x0,'CDF',[1 1],5); 
    x=iwt(a,d,an,dn,'CDF',[1 1]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(1,1) on Doppler noise')
   case 'demo2'   
    [a,d,an,dn]=wt(x0,'CDF',[1 3],5); 
    x=iwt(a,d,an,dn,'CDF',[1 3]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(1,3) on Doppler noise')
   case 'demo3'   
    [a,d,an,dn]=wt(x0,'CDF',[2 2],5); 
    x=iwt(a,d,an,dn,'CDF',[2 2]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(2,2) on Doppler noise')
   case 'demo4'   
    [a,d,an,dn]=wt(x0,'CDF',[2 4],5); 
    x=iwt(a,d,an,dn,'CDF',[2 4]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(2,4) on Doppler noise')
   case 'demo5'   
    [a,d,an,dn]=wt(x0,'CDF',[4 2],5); 
    x=iwt(a,d,an,dn,'CDF',[4 2]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(4,2) on Doppler noise')
   case 'demo6'   
    [a,d,an,dn]=wt(x0,'CDF',[6 8],5); 
    x=iwt(a,d,an,dn,'CDF',[6 8]);
    clf; plot(x{1}+x{2}+x{3}+x{4}+x{5}+x{6}-x0,'o'); ylim(minmax(x0))
    title('Reconstruction error for 5-level CDF(6,8) on Doppler noise')
   otherwise
    error('No such demo')
  end
end
