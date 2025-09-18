function [pte,sig,snr]=readOneDirection(map,Nra,Ndecl,ra,decl)

try decl; catch decl=ra.decl; ra=ra.ra; end;

skyPattern=SkyPattern(Nra,Ndecl);

dra=(mod(skyPattern(:,1)-ra+12,24)-12)/12*pi;
ddecl=(skyPattern(:,2)-decl)/180*pi;
metric=sqrt(dra.^2+ddecl.^2);

ind=find(metric==min(metric));

pte=mean(map(ind,1),1);
sig=mean(map(ind,2),1);
snr=pte./sig;
return
