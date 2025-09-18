function varargout=plotcont(c11,cmn,res,ofs)
% [axlim,handl,XYZ]=PLOTCONT(c11,cmn,res,ofs)
% 
% INPUT:
%
% c11   : [x,y]/[lon,lat] of upper left node
% cmn   : [x,y]/[lon,lat] of bottom right node
% res   : 0 regular coastlines (default) 
%       : 1 high-resolution coastlines and lakes
%       : 2 Global Mollweide projection
%       : 3 Three-dimensional coordinates
%       : 4 Three-dimensional coordinates, Northern hermisphere only
% ofs   : Longitude offset, e.g. 360 degrees
%
% Longitude ranges from 0 to 360; plots are centered on the PACIFIC.
%
% OUTPUT:
%
% axlim  : axis that fits snugly around the data.
% handl  : handle to the plotted line, and for Mollweide,
%          also the handle to the box around it
% XYZ    : the actual data points plotted (2D or 3D)
%
% AUSTRALIA:
% plotcont([90 10],[180 -60])
%
% WORLD:
% plotcont([0 90],[360 -90])
%
% Saved matrix as space-saving unsigned integer 
% - but that translates the NaN's into some  high number - take that out.
% Note how A==NaN does not return the NaN's!
%
% You'll have to maintain your own databases of continental outlines!
%
% Last modified by fjsimons-at-alum.mit.edu, 24.11.2004

%defval('c11',[90 10])
%defval('cmn',[180 -60])

defval('c11',[0 90])
defval('cmn',[360 -90])
defval('res',0)
defval('ofs',0)

pathname=getenv('IFILES');

switch res
 case 0
  fid=fopen([pathname '/COASTS/cont.mtl'],'r','b');
  cont=fread(fid,[5217 2],'uint16');
 case {1,3,4}
  fid=fopen([pathname '/COASTS/cost.mtl'],'r','b');  
  cont=fread(fid,[9598 2],'uint16');
 case 2
  load(fullfile(getenv('IFILES'),'GEOID','conm'))
  hold on
  handl{1}=plot(conxm,conym,'k');
  axlim=[];
  XYZ=[conxm conym];
  handl{2}=plot(xbox,ybox,'k');
end

if res~=2
  fclose(fid);
  cont=cont/100-90;
  cont(cont==max(max(cont)))=NaN;
  lon=cont(:,1); lat=cont(:,2);
  lon(~(lon>=c11(1) & lon<=cmn(1)))=NaN;
  lat(~(lat<=c11(2) & lat>=cmn(2)))=NaN;
  hold on
  handl=plot(lon+ofs,lat,'k','Linewidth',2);
  axlim=[min(lon([~isnan(lon) & ~isnan(lat)]))...
	 max(lon([~isnan(lon) & ~isnan(lat)]))...
	 min(lat([~isnan(lon) & ~isnan(lat)]))...
	 max(lat([~isnan(lon) & ~isnan(lat)]))];
  XYZ=[lon+ofs lat];
end
if res==3 | res==4
  lon=cont(:,1)/180*pi;
  lat=cont(:,2)/180*pi;
  rad=repmat(1.001,size(lat));
  [xx,yy,zz]=sph2cart(lon,lat,rad);
  XYZ=[xx yy zz];
  if res==4
    XYZ=XYZ(zz>0,:);
    % Now need to take out the annoying connecting lines
    % Distance between two consecutive points
    xx=XYZ(:,1); yy=XYZ(:,2); zz=XYZ(:,3);
    d=sqrt((xx(2:end)-xx(1:end-1)).^2+(yy(2:end)-yy(1:end-1)).^2);
    % D-level
    dlev=3;
    p=find(d>dlev*nanmean(d));
    % Now right after each of these positions need to insert a NaN;
    nx=insert(xx,NaN,p+1);
    ny=insert(yy,NaN,p+1);
    nz=insert(zz,NaN,p+1);
    XYZ=[nx(:) ny(:) nz(:)];
  end
end

vars={'axlim' 'handl' 'XYZ'};
for index=1:nargout
  varargout{index}=eval(vars{index});
end

hold off
