function [varargout]=SkyPattern(Nra,Ndecl,filename)

doSave=true;
try, filename; catch
  filename=['FullSkyPattern_',num2str(Nra),'_x_',num2str(Ndecl),'.txt'];
  doSave=false;
end;
rr=0:(Nra-1);
dd=0:(Ndecl-1);
data=zeros(Nra*Ndecl,2);
for r=rr
  for d=dd
    data(1+r*Ndecl+d,1)=r*24/Nra;
    data(1+r*Ndecl+d,2)=d*180/(Ndecl-1) - 90;
  end
end
if or(nargout == 0 , doSave)
  fprintf('Saving sky patterin in %s\n',filename);
  save(filename,'data','-ASCII');
end
if nargout == 1
  varargout{1}=data;
end
