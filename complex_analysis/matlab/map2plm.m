function [plm,lmcosi] = map2plm(map,RA,DECL,Lmax,cmplx)

% inverse of makemap.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try, cmplx;catch cmplx=0;end;

mapSize2=size(map,2);
if mapSize2 == 1

flagSkyPat=false;
if length(RA)>1
  Nra=length(unique(RA));
else
  Nra=RA;
  flagSkyPat=true;
end
if length(DECL)>1
  Ndecl=length(unique(DECL));
else
  Ndecl=DECL;
  flagSkyPat=true;
end    
if flagSkyPat
  pat  = SkyPattern(Nra,Ndecl);
  RA   = pat(:,1);
  DECL = pat(:,2);
end

  nRows=size(map,1);
  nPoints=Nra*Ndecl;
  if nRows == nPoints+1
    map=map(1:nPoints,:);
  elseif nRows ~= nPoints
    error('map has wrong dimensions');
    return;
  end


    %reshape the map
    ra   =reshape(RA    ,Ndecl,Nra);
    decl =reshape(DECL  ,Ndecl,Nra);
    Ymn  =reshape(map   ,Ndecl,Nra);

    % theta and phio vector
    theta = mean(ra,1)*15;
    phi   = mean(decl,2);

    % note that the phi goes from 90 to -90, while SkyPattern's default is from -90 to 90
    phi=transpose(flipud(phi));
    Ymn=flipud(Ymn);

    % normalize to 4 * pi, and add the 360deg line
    Ymn=Ymn*sqrt(4*pi);
    Ymn(:,end+1)=Ymn(:,1);
    theta(end+1)=360;
        
    [lmcosi,dw]=xyz2plm(real(Ymn),Lmax,'im');
    plm=plmreal2plm(lmcosi);
    
    if cmplx
      [lmcosi_i,dw]=xyz2plm(imag(Ymn),Lmax,'im');
      plm_i=plmreal2plm(lmcosi_i);
      plm=plm+i*plm_i;
    end   

else
  for kk=1:mapSize2
    [plm(:,kk),lmcosi] = map2plm(map(:,kk),RA,DECL,Lmax,cmplx);
  end
end



