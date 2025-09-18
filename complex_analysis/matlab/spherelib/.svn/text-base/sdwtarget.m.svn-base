function varargout=sdwtarget(co,TH,L,m,degres,SNm)
% [r,lon,lat,N,Nm,V,C]=SDWTARGET(co,TH,L,m,degres,SNm)
%
% Targets a region with a set of axisymmetric of Grunbaum's tapers
%
% INPUT
%
% co      [lon lat] coordinates of the center, in degrees
% TH      Radius of the window
% L       Bandwidth of the window
% m       Angular order of the window, -l=<m<=l
% degres  The degree resolution of the output field [default 1]
% SNm     The number of tapers you require [calculated default]
%
% OUTPUT:
%
% r       A spatial field with the taper right where you want it
%         in a matrix where the third dimension is the taper number
% lon     The longitudes at which r is evaluated, in degrees
% lat     The latitudes at which r is evaluated, in degrees
% N       The full Shannon number for this problem
% Nm      The partial Shannon number for this problem
% V       The calculated eigenvalues belonging to this set of tapers
% C       The spherical harmonic coefficients of the polar windows
%
% EXAMPLE:
% 
% sdwtarget('demo')
%
% Last modified by fjsimons-at-alum.mit.edu, 28.3.2005

if ~isstr(co) 
  
  defval('co',[0 0])
  defval('TH',20)
  defval('L',18)
  defval('m',0)
  defval('degres',1)

  if(m>L)
    error('Order cannot exceed degree')
  end
  
  % Calculate the taper coefficients centered on the pole
  [E,Vg,th,C,T,V]=grunbaum(TH,L,m,0);
  
  % Calculate the full Shannon number
  N=(L+1)^2*(1-cos(TH/180*pi))/2;

  % Calculate the asymptotic partial Shannon number
  [Nm,Nsum]=nsubm(N,L,1,L);
  
  % Supply a default number of required tapers
  defval('SNm',round(Nm(abs(m)+1)))
  
  % This is for the output only
  Nm=Nm(abs(m)+1);

  % Make a blank array of indices and coefficients
  [em,el,mzero]=addmon(L);
  Cb=[el em repmat(0,length(em),2)];

  % Make a blank window array 
  r=repmat(NaN,[180/degres+1 360/degres+1 SNm]);

  % Find the appropriate rotation angles
  alpha=0; % Around z axis
  beta=90-co(2); % To a desired colatitude
  gamma=co(1); % To a desired longitude

  % Do the following for every taper
  for index=1:SNm
    % Put the zonal coefficients in the right place
    Cb(mzero(abs(m)+1:end)+abs(m),3+(m<0))=C(:,index);
    
    % Rotate the tapers to a position on the sphere
    Cr=plm2rot(Cb,alpha,beta,gamma);

    % Calculate the spatial functions
    [r(:,:,index),lon,lat]=plm2xyz(Cr,degres);
    
    % Adjust the sign
    if max(max(abs(r(:,:,index))))==abs(min(min(r(:,:,index))))
      r(:,:,index)=-r(:,:,index);
    end
  end

  % Output
  varn={'r','lon','lat','N','Nm','V','C'};
  for index=1:nargout
    varargout{index}=eval(varn{index});
  end
else
  col=round(rand*180);
  lon=round(rand*360); 
  r=sdwtarget([lon 90-col],20,12,0,1,4);
  clf
  for index=1:4
    ah(index)=subplot(2,2,index);
    imagef([],[],r(:,:,index))
    set(gca,'xtick',[0 lon 360],'xtickl',[0 lon 360])
    set(gca,'ytick',[-90 90-col 90],'ytickl',[-90 90-col 90])
  end
  supertit(ah(1:2),sprintf('Lon= %i; Lat= %i',lon,90-col));
  movev(ah(1:2),-0.025)
  longticks(ah)
  set(ah,'xgrid','on','ygrid','on')
end
