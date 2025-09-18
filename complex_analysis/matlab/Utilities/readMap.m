function [map,coord1,coord2]=readMap(filename,fileType,numBins,InjectAsSpH,conversionLmax,conversionDeg)

%
%  readMap --- load a map (spherical harmonics or pixel) for injection
%
%  parameters:
%    filename      : filename
%    fileType      : 0: pixel Map; 1: real SpH; 2: complex SpH
%    numBins       : number of frequency bins in file
%                    Optional - default 1
%    InjectAsSpH   : output is spherical harmonics if true
%                              pixel map if false
%                    Optional - default true
%  only required when is conversion is done
%    conversionLmax: for pixel -> SpH: Lmax
%    conversionDeg : for SpH -> pixel: resolution in degree 
%                   
%
%    output:         depending on InjectAsSpH it can be either of two
%                    InjectAsSpH=true:
%    map           : vector of complex spherical harmonics
%    coord1        : vector of l
%    coord2        : vector of m
%                    InjectAsSpH=false:
%    map           : vector of pixel values
%    coord1        : vector of RA
%    coord2        : vector of DECL
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@caltech.edu
%
%  $Id: readMap.m,v 1.1 2008-01-30 19:29:49 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try, numBins; catch numBins=1; end;
try, InjectAsSpH; catch InjectAsSpH=true; end;

data=load(filename);
sz=size(data);

switch fileType
  case 0,
    hasRADEC = checkLength(sz(2),numBins);
    if ~hasRADEC
      error('Pixel map needs RA and DECL column.');
    end
    if InjectAsSpH
      for ii=1:numBins
        map(:,ii) = map2plm(data(:,2+ii),data(:,1),data(:,2),conversionLmax);
      end
      [coord1,coord2]=getLMvec(conversionLmax);
    else
      coord1 =data(:,1);
      coord2 =data(:,2);
      map    =data(:,3:end);
    end
    
  case 1,
    hasLM = checkLength(sz(2),2*numBins);
    lmax=(sqrt(8*sz(1)+1)-3)/2;
    if mod(lmax,1)~=0
      error('Number of rows not consistent with real spherical harmonics');
    end
    if numBins==1
      [plm,coord1,coord2] = plmreal2plm(data);
    else
      for ii=1:numBins
        if hasLM
          [plm(:,ii),coord1,coord2] = plmreal2plm([data(:,1:2),data(:,(1:2)+2*ii)]);
	else
          [plm(:,ii),coord1,coord2] = plmreal2plm(             data(:,(1:2)+2*(ii-1)));
	end
      end
    end  
    if InjectAsSpH
      map=plm;      
    else
      for ii=1:numBins
        [map(:,ii),coord1,coord2,dOmega]=makemap(plm(:,ii),conversionDeg);
      end
    end

  case 2,
    hasLM = checkLength(sz(2),2*numBins);
    lmax=sqrt(sz(1))-1;
    if mod(lmax,1)~=0
        error('Number of rows not consistent with complex spherical harmonics');
    end
    if hasLM
      plm=data(:,3:2:end)+i*data(:,4:2:end);
      coord1=data(:,1);
      coord2=data(:,2);
    else
      plm=data(:,1:2:end)+i*data(:,2:2:end);
      [coord1,coord2]=getLMvec(lmax);
    end
    if InjectAsSpH
      map=plm;      
    else
      for ii=1:numBins
        [map(:,ii),coord1,coord2,dOmega]=makemap(plm(:,ii),conversionDeg);
      end
    end

  otherwise,
    error('File type unknown: 0: pixel Map; 1: real SpH; 2: complex SpH');
end

return;

function hasLM = checkLength(l,N)
  if l==N
    hasLM=false;
  elseif l==N+2
    hasLM=true;
  else
     error('Wrong number of columns for loading as spherical harmonics.');
  end
return




