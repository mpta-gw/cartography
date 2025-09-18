function writeMap(plm,filename,fileType,deg)

%
%  readMap --- wrtie a map (spherical harmonics) to file
%
%  parameters:
%    plm        : complex spherical harmonics
%    filename   : filename
%    fileType   : 0: pixel Map; 1: real SpH; 2: complex SpH
%    deg        : resolution of pixel map in degree
%                   only for fileType=0; default = 1 
%
%   ordering in plm:
%              |  plm
%   ---------------------------------
%   M=-L,L=Lmax|
%   .......    |
%   M=-1,L=Lmax|
%   .......    |
%   M=-1,L=1   |
%   M=0,L=0    |
%   M=0,L=1    |
%   .......    |
%   M=0,L=Lmax |
%   M=1,L=1    |
%   .......    |
%   M=1,L=Lmax |
%   .......    |
%   M=L=Lmax   |
%
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@caltech.edu
%
%  $Id: writeMap.m,v 1.1 2008-01-30 19:29:49 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz=size(plm);
lmax=sqrt(sz(1))-1;

switch fileType
  case 0,
    try, deg; catch deg=1; end;
    for ii=1:sz(2)
      [map(:,ii),RA,DECL,dOmg] = makemap(plm(:,ii),deg);
    end
    data=[RA,DECL,map];
  case 1,
    data=plm2plmreal(plm(:,1));
    if sz(2)>1
      for ii=2:sz(2)
        plmreal=plm2plmreal(plm(:,ii));
	data(:,(1:2)+2*ii)=plmreal(:,3:4);
      end
    end
  case 2,
    [lvec,mvec]=getLMvec(lmax);
    data(:,1:2)=[lvec,mvec];
    for ii=1:sz(2)
      data(:,(1:2)+2*ii)=[real(plm(:,ii)),imag(plm(:,ii))];
    end
  otherwise,
    error('File type unknown: 0: pixel Map; 1: real SpH; 2: complex SpH');
end

save(filename,'data','-ASCII','-DOUBLE');

