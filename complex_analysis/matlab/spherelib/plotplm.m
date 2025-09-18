function varargout=plotplm(data,lon,lat,meth,degres,th0)
% PLOTPLM(data,lon,lat,meth,degres,th0)
% PLOTPLM(lmcosi,[],[],meth,degres,th0)
% PLOTPLM(l,m,cosi,meth,degres,th0)
% [data,ch,ph]=PLOTPLM
%
% Global plotting routine for fields coming out of PLM2XYZ,
% plotting REAL spherical harmonics.
%
% INPUT:
%
% data           2-D data/or LMCOSI matrix
% lmcosi         Degree, Order, Cosine, Sine Coefficient matrix
%                (column size is either 4 or 6)
% l              Angular degree
% m              Angular order
% cosi           Sine or cosine component
% lon,lat        2-D grid in radians [default (0,2pi) and (-pi/2,pi/2)]
% meth           1  Mollweide projection
%                2  3-D sphere, no topography
%                3  Plots the spectrum
%                4  Flat rectangular projection
%                5  Plot looking down on the North Pole
% degres         Resolution in degree (checked for Nyquist)
% th0            Small circle at colatitude th0, in degrees
%
% 'lon' and 'lat' are in radians.
%
% OUTPUT:
%
% data           Spatial data in case coefficients were given
%                Handles to spectral log-log plots for meth 3
% ch             Handle to the continents
% ph             Handle to the plates or the small circles
%
% See also PLOTONSPHERE, PLOTONEARTH, CPX2RSH, RSH2CPX.
%
% Last modified by fjsimons-at-alum.mit.edu, 25.05.2005

defval('meth',1)
defval('degres',1)
defval('th0',[])

if (size(data,2)==4 | size(data,2)==6) & meth<3
  [data,lon,lat]=plm2xyz(data,degres);
  lon=lon*pi/180;
  lat=lat*pi/180;
elseif length(data)==1
  defval('lat',1)
  [m,l,mz,lmcosi]=addmon(data);
  lmcosi(mz(data+1)+lon,2+lat)=1;
  [data,lon,lat]=plm2xyz(lmcosi,degres);
  lon=lon*pi/180;
  lat=lat*pi/180;
end

switch meth
 case 1
   defval('lon',linspace(0,2*pi,size(data,2)))
   defval('lat',linspace(pi/2,-pi/2,size(data,1)))
   % Need to make special provisions for PCOLOR compared to IMAGESC
   % The input values are PIXEL centered
   dlat=(lat(1)-lat(2))/2;
   dlon=(lon(2)-lon(1))/2;
   lat=[lat(1) lat(2:end)+dlat lat(end)];
   lon=[lon(1) lon(2:end)-dlon lon(end)];
   [lon,lat]=meshgrid(lon,lat);
   [xgr,ygr]=mollweide(lon,lat,pi);
   pc=pcolor(xgr,ygr,adrc(data)); shading flat
   [ax1,p1,XY1]=plotcont([],[],2);
   [ax2,XY2]=plotplates([],[],2);
   axis image
   axis off
 case 2
  % Make sphere for the data
  defval('lon',linspace(0,360,100)/180*pi);
  defval('lat',linspace(90,-90,100)/180*pi);
  [lon,lat]=meshgrid(lon,lat);
  rads=ones(size(lat));
  [x,y,z]=sph2cart(lon,lat,rads);
  surface(x,y,z,'FaceColor','texture','Cdata',data);   
  hold on
  % Plot the continents
  [axlim,handl,XYZ]=plotcont([0 90],[360 -90],3);  
  delete(handl)
  hold on
  p1=plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'k-','LineWidth',1.5);
  axis image
  shading flat
  view(140,30)
  set(gca,'color',[.9 .9 .9]-.1)
  axis('off')
  colormap default
  ax2=[];
 case 3
  [sdl,l,bta,lfit,logy,logpm]=plm2spec(data,2);
  % The following is now well taken care of by LIBBRECHT
  % l=l(l<256 & l>0); 
  % sdl=sdl(l<256 & l>0);
  l=l(l>0);
  sdl=sdl(l>0);
  disp(sprintf('PLOTPLM Spectrum plot from %i to %i',...
	       l(1),l(end)))
  clear data
  data(:,1)=loglog(l,sdl);
  set(data(:,1),'Marker','+')
  grid on
  tix=sort([l(1) 10.^[1 2 3] l(end)]);
  tix=tix(tix<=l(end)); 
  xlim([l(1) l(end)])
  set(gca,'xtick',tix)
  tiy=round(minmax(log10(sdl(~~sdl))));
  set(gca,'ytick',10.^[tiy(1):1:tiy(2)],...
	  'yminortick','off','yminorgrid','off')
  hold on
  data(:,2)=loglog(lfit,logy);
  % Do not use axis tight or openup on loglog
  ylim([0.9 1.1].*minmax([sdl(:) ; logy(:)]))  
  set(data(:,2),'Color','r')
  longticks(gca)
  data(:,3)=title(sprintf('%s= %8.3f','\beta',bta));
  data(:,4)=bta;
 case 4
   imagef([0 90],[360 -90],data)
   [ax1,p1,XY1]=plotcont([],[],1);
   [ax2,XY2]=plotplates([],[],1);
   axis image
   axis off  
 case 5
  % Project northern hemisphere only
  lon=linspace(0,2*pi,size(data,2));
  lat=linspace(pi/2,0,floor(size(data,1)/2));
  [LON,LAT]=meshgrid(lon,lat);
  % Radius from 0 to 1; longitude is azimuth
  r=cos(LAT);
  x=r.*cos(LON);
  y=r.*sin(LON);
  % Only project the upper hemisphere
  X=linspace(-1,1,500);
  Y=linspace(-1,1,500);
  warning off MATLAB:griddata:DuplicateDataPoints
  % Watch out: the POLE is a new problem with the latest version of
  % Matlab; need to fake this entirely
  y(1,:)=y(2,:)/2;
  x(1,:)=x(2,:)/2;
  Z=griddata(x,y,data(1:length(lat),:),X,Y(:));
  disp('Data gridded by PLOTPLM')
  warning on MATLAB:griddata:DuplicateDataPoints
  % imagefnan([-1 1],[1 -1],Z,gray(10),minmax(Z(:)))  
  % colormap(gray(10)); hold on
  imagefnan([-1 1],[1 -1],Z,kelicol,minmax(Z(:)))
  colormap(kelicol); hold on
  ax2(1)=circ(1); ax2(2)=circ(sin(th0*pi/180));
  set(ax2(1:2),'LineW',1)
  set(ax2(2),'LineS','--')
  axis([-1.0100    1.0100   -1.0100    1.0100])
  axis off
  [axlim,handl,XYZ]=plotcont([0 90],[360 -90],4);
  delete(handl)
  hold on
  p1=plot(XYZ(:,1),XYZ(:,2),'k-','LineWidth',1);
  data=Z;
 otherwise
  error('Not a valid method')
end

varn={'data','p1','ax2'};
for index=1:nargout
  varargout{index}=eval(varn{index});
end


