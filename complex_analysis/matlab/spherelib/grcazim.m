function [lon2,lat2]=grcazim(lonlat,delta,angl)
% [lon2,lat2]=GRCAZIM([lon1 lat1],delta,angl)
%
% Calculated location of points lying on a great circle through a certain
% point with a certain azimuth at a certain epicentral distance. One of
% the aguments can be a vector, OR both [lon lat] and delta can be
% vectors of the same length.
%
% INPUT:
%
% [lon1 lat1]        First point (can be vector), in degrees
% delta              Epicentral distance, in km (can be vector)
% angl               Angle (vector) WEST from NORTH, in degrees
%                    Note: this is true course rather than azimuth
%                    NORTH-SOUTH: angl=0,  and delta>0 is NORTH
%                    EAST-WEST:   angl=90, and delta>0 is WEST
%
% OUPUT:
%
% [lon2 lat2]        Coordinates of the great circle
%
% EXAMPLE:
%
% plotcont([90 10],[180 -60])
% hold on
% [lon2,lat2]=grcazim([134 -24],[-1000:5:1000],rand*90); 
% plot(134,-24,'bs'); plot(lon2,lat2,'r')
% [lon2,lat2]=grcazim([134 -24],1000,linspace(0,360,20)); 
% plot(lon2,lat2,'b')
% [lon2,lat2]=grcazim([124 -20 ; 114 -30],1000,linspace(0,360,20)); 
% plot(lon2,lat2,'b')
% [lon2,lat2]=grcazim([124 -20 ; 114 -30],0:100:1000,45); 
% plot(lon2,lat2,'g')
%
% Last modified by fjsimons-at-alum.mit.edu, Mar 7th, 2004

lon1=lonlat(:,1); delta=delta(:);
lat1=lonlat(:,2); angl=angl(:);
lon1=lon1*pi/180; lat1=lat1*pi/180;
delta=delta*1000/fralmanac('Radius');
angl=angl*pi/180;

% Only one vector at a time except if lonlat and delta have the same length
if (length(lat1)==1 & length(delta)==1) | (length(lat1)==1 & length(angl)==1)
  lat2 =asin(sin(lat1)*cos(delta)+cos(lat1)*sin(delta)*cos(angl));
  dlon=atan2(sin(angl)*sin(delta)*cos(lat1),cos(delta)-sin(lat1)*sin(lat2));
  lon2=mod(lon1-dlon+pi,2*pi)-pi;
elseif length(delta)==1
  lat2=asin(repmat(sin(lat1)*cos(delta),1,length(angl))+...
	     cos(lat1)*sin(delta)*cos(angl(:)'));
  dlon=atan2(cos(lat1)*sin(angl(:)')*sin(delta),...
	     cos(delta)-repmat(sin(lat1),1,length(angl)).*sin(lat2));
  lon2=mod(repmat(lon1,1,length(angl))-dlon+pi,2*pi)-pi;
  lat2=lat2';
  lon2=lon2';
elseif length(angl)==1
  lat2=asin(sin(lat1)*cos(delta(:)')+cos(lat1)*sin(delta(:)')*cos(angl));
  dlon=atan2(cos(lat1)*sin(angl)*sin(delta(:)'),...
	     repmat(cos(delta(:)'),length(lat1),1)-...
	     repmat(sin(lat1),1,length(delta)).*sin(lat2));
  lon2=mod(repmat(lon1,1,length(delta))-dlon+pi,2*pi)-pi;
  lat2=lat2';
  lon2=lon2';
elseif length(lat1)==length(delta)
  for index=1:length(lat1)
    lat2(:,index) =asin(sin(lat1(index))*cos(delta(index))+...
			cos(lat1(index))*sin(delta(index))*cos(angl));
    dlon=atan2(sin(angl)*sin(delta(index))*cos(lat1(index)),...
	       cos(delta(index))-sin(lat1(index))*sin(lat2(:,index)));
    lon2(:,index)=mod(lon1(index)-dlon+pi,2*pi )-pi;
  end
end

i=lon2<0;

lon2(i)=lon2(i)+2*pi;

lon2=lon2*180/pi; 
lat2=lat2*180/pi;


