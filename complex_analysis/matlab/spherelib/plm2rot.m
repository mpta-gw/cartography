function varargout=plm2rot(lmcosi,alp,bta,gam,method,rlcp)
% [lmcosip,spec1,spec2]=PLM2ROT(lmcosi,alp,bta,gam,method,rlcp)
%
% Rotates spherical harmonics by rearranging their
% coefficients after computing the rotation matrix.
%
% INPUT:
%
% lmcosi           Matrix with [l m Ccos Csin] in order, m>=0
% alp, bta, gam    Three Euler angles, in degrees
%                  alp (0<2pi): rotates around z increasing from x to y
%                  bta (0<pi):  rotates around new y increasing from z to x
%                  gam (0<2pi): rotates around old z increasing from x to y
% method           'dlmb' using DLMB [default]
%                  'blanco' using BLANCO with bta=90
% rlcp             1 if coefficients belong to REAL harmonics [default]
%                  2 if coefficients belong to COMPLEX harmonics 
%
% OUTPUT:
%
% lmcosip          Matrix with [l m Ccosp Csinp] rotated
% spec1            Spectral density of input field
% spec2            Spectral density of rotated field
%
% EXAMPLES:
%
% plm2rot('demo1') % For a random input with some spectral slope
% plm2rot('demo2') % For Master's example coefficient set
% plm2rot('demo3') % For Wieczorek-Simons taper 
% plm2rot('demo4') % For a random pure spherical harmonic
% plm2rot('demo5') % For a random non-azimuthal cap
% plm2rot('demo6') % Difference between BLANCO and DLMB
% plm2rot('demo7') % Yet another example
% plm2rot('demo8') % Presentation-quality figure
% 
% Last modified by fjsimons-at-alum.mit.edu, 23.2.2005

if ~isstr(lmcosi)
  defval('alp',0)
  defval('bta',0)
  defval('gam',0)
  defval('method','dlmb')
  defval('rlcp',1)
  
  if alp==0 & bta==0 & gam==0
    lmcosip=lmcosi;
    spec1=[]; spec2=[];
    varn={'lmcosip','spec1','spec2'};
    for index=1:length(nargout)
      varargout{index}=eval(varn{index});
    end
    return
  end
  
  disp(sprintf('Rotation over %i %i %i',...
	       round(alp),round(bta),round(gam)))
  
  alp=alp*pi/180-pi/2;
  bta=-bta*pi/180;
  gam=gam*pi/180+pi/2;
  
  % Determine maximum L  
  L=max(lmcosi(:,1));
  % Initialize rotated array 
  lmcosip=lmcosi;

  % Methods to calculate the beta-rotation matrix
  switch method
   case 'dlmb'
    % Use Master's method
    disp('Using DLMB method')
    D=dlmb(L);
   case 'blanco'
    disp('Using Blanco method')
    warning('Not recommended')
    % Use algorithm of Blanco, Florez and Bermejo (1997).
    D=blanco(L,90);
   case 'blanco2'
    [D,DC]=blanco(L,bta*180/pi);
   otherwise
    error('Specify valid method')
  end

  spec1=plm2spec(lmcosi);

  if rlcp==1
    lmcosi=rsh2cpx(lmcosi);
  end

if strmatch('blanco2',method)
  error('Not supported')
  for l=0:L
    [C,b,e]=shcos(lmcosi,l);
    S=shsin(lmcosi,l);
    cpl=[C ; S(2:end)];
    cplp=DC{l+1}*cpl;
    lmcosip(b:e,3)=cplp(1:l+1);
    lmcosip(b:e,4)=[0 ; cplp(l+2:end)];
  end
  if rlcp==1
    lmcosip=cpx2rsh(lmcosip);
  end
  varargout{1}=lmcosip;
  spec2=plm2spec(lmcosip);
  dif=abs(sum(spec1-spec2))/L;
  disp(sprintf('Average spectral error %8.3e',dif))
  if dif>1e-10
    error('Something wrong - real/complex mixup?')
  end
  return
end
  % ALPHA rotation
  [C,S]=rotcof(lmcosi(:,3),lmcosi(:,4),alp);
  
  % Loop over all degrees l
  % BETA rotation
  cbet=cos([0:L]*bta)';
  sbet=sin([0:L]*bta)';
  for l=0:L
    li=l+1;
    % Loop over all orders m
    % Make a custom alternating matrix with zeros on the diagonal
    [i,j]=meshgrid(1:li,1:li);
    IC=mod((i+j)+mod(li,2),2)*2;  IC(:,1)=1;
    IS=~mod((i+j)+mod(li,2),2)*2; IS(:,1)=1;

    Cp=D{li}'.*IC*shcos(C,l);
    Sp=D{li}'.*IS*shsin(S,l);

    Cpp=Cp.*cbet(1:l+1)-Sp.*sbet(1:l+1);
    Spp=Sp.*cbet(1:l+1)+Cp.*sbet(1:l+1); 

    Cpp=D{li}.*IC*Cpp;
    Spp=D{li}.*IS*Spp;

    % Fill up the matrix
    b=addmup(l-1)+1;
    e=addmup(l);

    Crot(b:e,1)=Cpp;
    Srot(b:e,1)=Spp;
  end

  % GAMMA Rotation
  [Cp,Sp]=rotcof(Crot,Srot,gam);
  lmcosip(:,3)=Cp;
  lmcosip(:,4)=Sp;
  if rlcp==1
    lmcosip=cpx2rsh(lmcosip);
  end
  % Verify that spectrum remains unchanged
  spec2=plm2spec(lmcosip);
  dif=abs(sum(spec1-spec2))/L;
  disp(sprintf('Average spectral error %8.3e',dif))
  if dif>1e-10
    error('Something wrong - real/complex mixup?')
  end
  % Output
  varn={'lmcosip','spec1','spec2'};
  for index=1:length(nargout)
    varargout{index}=eval(varn{index});
  end
else
  demonr=lmcosi;
  switch demonr
   case 'demo1'
    lmax=ceil(20);
    [m,l,mzero]=addmon(lmax);               
    c=randn(addmup(lmax),2).*([l l].^(-1)); 
    c(mzero,2)=0; c(1,1)=1;
    lmcosi=[l m c]; 
    th0=NaN;
   case 'demo2'
    lmcm=load(fullfile(getenv('IFILES'),'masters')); 
    lmcosi=lmcm.mas; l=6; % These are complex
    lmcosi=cpx2rsh(lmcosi);
    th0=NaN;
   case 'demo3'
    th0=40;
    [E,V,Lmax,TH,C]=wieczorek(th0); 
    [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)];
    lmcosi(mzero,3)=C(:,1);
    th0=NaN;
   case 'demo4'
    Lmax=30;
    [em,el,mz]=addmon(Lmax); M=round(rand*Lmax);
    c=zeros(length(el),2);
    c(mz(end)+M(1),1)=1;
    lmcosi=[el em c];
    th0=NaN;
   case 'demo5'
    th0=40;
    [E,V,Lmax,TH,C]=wieczorek(th0); [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)]; 
    lmcosi(mzero,3)=C(:,1); 
%    lmcosi(mzero+1,3)=C(:,1);
    lmcosi(mzero,4)=0;
   case 'demo6'
    lmax=ceil(20);
    [m,l,mzero]=addmon(lmax);               
    c=randn(addmup(lmax),2).*([l l].^(-1)); 
    c(mzero,2)=0; c(1,1)=1;
    lmcosi=[l m c]; 
    alp=round(rand*360);
    bta=round(rand*180); 
    gam=round(rand*360); 
    lmcosip1=plm2rot(lmcosi,alp,bta,gam,'dlmb');
    lmcosip2=plm2rot(lmcosi,alp,bta,gam,'blanco');
    dif=abs(sum(sum(lmcosip1(:,3:4)-lmcosip2(:,3:4))));
    disp(sprintf('BLANCO vs DLMB %8.3e',dif))
    if dif>1e-10
      error('Something wrong - correct beta angle?')
    end
    return
   case 'demo7'
    th0=40;
    [E,V,Lmax,TH,C]=wieczorek(th0); [m,l,mzero]=addmon(Lmax);
    lmcosi=[l m zeros(length(l),2)];
    lmcosi(mzero+1,4)=C(:,1);
    lmcosi(mzero,4)=0;
   case 'demo8'
    clf
    th0=40;
    [E,V,Lmax,TH,C]=wieczorek(th0); 
    [m,l,mzero]=addmon(Lmax);
    alp=30;
    bta=110;
    gam=270;
    ah=krijetem(subnum(2,2));
    lmcosi=[l m zeros(length(l),2)];
    [xgb,ygb]=mollweide(gam*pi/180,(90-bta)*pi/180,pi);
    [xa,ya]=longitude(gam,3);
    [xo,yo]=latitude(90-bta,3);
    [xi,yi]=caploc([gam 90-bta],th0,[],2);
    for index=1:4
      lmcosi(mzero,3)=C(:,index);
      lmcosip=plm2rot(lmcosi,alp,bta,gam,'dlmb');
      th0=NaN;
      axes(ah(index))
      [ch(index),ph(index)]=plott3(lmcosip);
      hold on
      ps(index)=plot(xgb,ygb,'s');
      ha(index)=plot(xa,ya,'r');
      ho(index)=plot(xo,yo,'r');
      hi(index)=plot(xi,yi,'r');            
    end   
    set(ah,'CameraV',6)
    movev(ah(1:2),-.15)
    keyboard
    return
   otherwise
    error('Not a valid demo.')
  end
  clf
  % Now perform random rotation
  alp=round(rand*360);
  bta=round(rand*180);
  gam=round(rand*360);
  lmcosip1=plm2rot(lmcosi,alp,bta,gam,'dlmb');

  % And plot the shit up
  ah=krijetem(subnum(4,2));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(1))
  [data,labs(1,:)]=plott1(lmcosi);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(2))
  [data2,labs(2,:)]=plott1(lmcosip1);
  seemax(ah(1:2),[1 2 4])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(3))
  [h(1),labs(3,:)]=plott2(data);
  [x,y,z]=latitude(90-th0,1);
  ha(2)=plot3(x,y,z,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axes(ah(4))
  [h(2),labs(4,:)]=plott2(data2);
  hold on
  [x,y,z]=sph2cart(gam/180*pi,(90-bta)/180*pi,1.01);
  ps(1)=plot3(x,y,z,'s');
  [x,y,z]=latitude(90-bta,1);
  h(16)=plot3(x,y,z,'w-');
  [x,y,z]=longitude(gam,1);
  h(11)=plot3(x,y,z,'w-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(5))
  [h(3),h(4)]=plott3(data);
  hold on
  [x,y,z]=latitude(90-th0,3);
  ha(3)=plot(x,y,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(6))
  [h(5),h(6)]=plott3(data2);
  hold on
  [xgb,ygb]=mollweide(gam*pi/180,(90-bta)*pi/180,pi);
  ps(4)=plot(xgb,ygb,'s');
  [x,y]=longitude(gam,3);
  h(12)=plot(x,y,'w');
  [x,y]=latitude(90-bta,3);
  h(13)=plot(x,y,'w');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(7))
  [h(7),h(8)]=plott4(data);
  hold on
  [x,y,z]=latitude(90-th0,2);
  ha(1)=plot(x,y,'y-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  axes(ah(8))
  [h(9),h(10)]=plott4(data2);
  hold on
  ps(5)=plot(gam,90-bta,'s');
  [x,y,z]=latitude(90-bta,2);
  h(14)=plot(x,y,'w-');
  [x,y,z]=longitude(gam,2);
  h(15)=plot(x,y,'w-');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  longticks(ah(7:8))
  set(ah(7:8),'xtick',[0:90:360],'ytick',[-90:45:90])
  deggies(ah(7:8))
  set(ah(1:2),'camerav',8.5)
  set(ah(3:4),'camerav',8.5)
  set(ah(5:6),'camerav',5)
  set(h,'LineW',0.5)
  set(ps([1 4 5]),'MarkerF','g','MarkerE','k')
  t(1)=supertit(ah(1:2),sprintf('%s= %i ; %s= %i ; %s= %i',...
				'\alpha',alp,'\beta',bta,'\gamma',gam));
  fig2print(gcf,'tall')
  figdisp([],demonr)
  set(gcf,'inv','off','color','w')
end

% Auxiliary plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spat,xyz]=plott1(cofs)
spat=plotplm(cofs,[],[],2); axis image; hold off
plotonsphere(spat,0.2); axis image on; 
xyz(1)=xlabel('x');
xyz(2)=ylabel('y');
xyz(3)=zlabel('z'); grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,xyz]=plott2(spat)
ch=plotonearth(spat,1); axis image on; box on ; grid on
set(ch,'Color','w')
xyz(1)=xlabel('x');
xyz(2)=ylabel('y');
xyz(3)=zlabel('z'); grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,ph]=plott3(spat)
[jk,chh,ph]=plotplm(spat,[],[],1); axis image; hold off
ch=chh{1}; set([ch ph],'Color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,ph]=plott4(spat)
imagef([0 90],[360 -90],spat)
[jk,ch]=plotcont;
ph=plotplates;
set([ch ph],'Color','w')
