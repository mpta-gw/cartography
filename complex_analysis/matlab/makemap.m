function [map,RA,DECL,dOmg] = makemap(plm,deg,norm,cmplx)

%  Input:
%           plm  -  complex spherical harmonics packed as:
%
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
%           deg  -  desired resolution of map
%           norm -  1=> power normalized to dOmega (for summing over point
%           sources)  2=> power at each point (for input to plotMapAitoff,
%           for comparing with the output of plmToPlmPix.m)
%
%  Output:
%           map  -  a pixelated map of the spherical harmonic.  First
%           column: right ascension in hours [0, 24).  Second column: declination
%           in degrees [-90, 90].  Third column: power*dOmega (if norm=1), 
%           or power (if norm=0) such that
%           (the sum over all the points of (power^2 * dOmega)) = 1
%
%
%  Routine written by Madeleine Udell.
%  Some corrections by Stefan Ballmer
%  Contact madeleine_udell@yale.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try, norm; catch norm=0; end;
try, deg;  catch deg =1; end;
try, cmplx;catch cmplx=0;end;

plmSize2=size(plm,2);
if plmSize2 == 1
%  convert to real spherical harmonics
    lmcosi=plm2plmreal(plm);

%  shift the locations of the declination points so that none are at the
%  poles for higher accuracy with a given resolution

    % first make a map with twice the resolution
    [Ymn, theta, phi] = plm2xyz(lmcosi, deg);
    
    if cmplx
      disp("Here")
      lmcosi_i=plm2plmreal(-i*plm);
      Ymn_i = plm2xyz(lmcosi_i, deg);
      Ymn=Ymn+i*Ymn_i;
    end
         
    % plm2xyz gives theta and phi in degres, phi (declination) from -90 to 90 
    Nra=length(theta) - 1;
    Ndecl=length(phi);
    map=zeros(Nra*Ndecl,3);
    
    %chop the 360deg line, and normalize to 1
    Ymn=Ymn(:,1:end-1) ./sqrt(4*pi);
    theta=theta(1:end-1);

    % note that the phi goes from 90 to -90, while SkyPattern's default is from -90 to 90
    phi=flipud(transpose(phi));
    Ymn=flipud(Ymn);
           
    % define ra and decl matrices
    ra   = ones(Ndecl,1)*theta/15;
    decl = phi*ones(1,Nra);

    % does the map already contain the factor dOmega ? Scale it out!
    dOmega=cos(decl/180*pi)* ( 2*pi^2/(Nra*(Ndecl-1)) );
    if norm
      Ymn=Ymn.*dOmega;
    end
    
    %reshape the map
    RA   =reshape(ra    ,Nra*Ndecl,1);
    DECL =reshape(decl  ,Nra*Ndecl,1);
    map  =reshape(Ymn   ,Nra*Ndecl,1);
    dOmg =reshape(dOmega,Nra*Ndecl,1);
    
    % calculate norm
    if norm
      nor=sum(map'*(map./dOmg));
    else
      nor=sum(map'*(map.*dOmg));
    end
    %disp(['Norm of map: ',num2str(nor)]);
else
  for kk=1:plmSize2
    [map(:,kk),RA,DECL,dOmg] = makemap(plm(:,kk),deg,norm,cmplx);
  end
end



