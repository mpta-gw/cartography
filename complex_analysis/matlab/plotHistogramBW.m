function plotHistogramBW(map,Nra,Ndecl,Nindep,binWidth)
defval('Nindep',100);
defval('binWidth',0.5);

  nRows=size(map,1);
  nPoints=Nra*Ndecl;
  if nRows == nPoints+1
    map=map(1:nPoints,:);
  elseif nRows ~= nPoints
    error('map has wrong dimensions');
    return;
  end

set(0,'DefaultLineLineWidth', 2.0);
pat=SkyPattern(Nra,Ndecl);

dphi=360/Nra;
dtheta=180/(Ndecl-1);
pixelarea=dtheta*dphi*cos(pat(:,2)/180*pi);

snr=map(:,1)./map(:,2);

numLoc=length(snr);

bmin=min(-5,min(snr));
bmax=max(+5,max(snr));
bins=bmin:binWidth:bmax;
fine=bmin:binWidth/10:bmax;
N=hist([],bins);
for d=pat(1:Ndecl,2)'
  ind=find(pat(:,2)==d);
  M=hist(snr(ind),bins);
  N=N+M*dtheta*dphi*cos(d/180*pi);
end;
figure;
bar(bins,N);
set(gca,'FontSize',16);
set(gca,'FontWeight','bold');
set(gca,'LineWidth', 2.0);
set(gca,'GridLineStyle','--');
set(gca,'XColor',[0 0 0]);
set(gca,'YColor',[0 0 0]);

y=1./sqrt(2*pi).*exp(-bins.^2/2)*(2/pi*numLoc)*(bins(2)-bins(1));

x0=sum(N.*bins)/sum(N);
sigma=sqrt(sum(N.*(bins-x0).^2)/sum(N));
z=1./sqrt(2*pi*sigma^2).*exp(-(fine-x0).^2/(2*sigma^2))*(2/pi*numLoc)*(bins(2)-bins(1));
zideal=1./sqrt(2*pi*1^2).*exp(-(fine-0).^2/(2*1^2))*(2/pi*numLoc)*(bins(2)-bins(1));

hold on;
leg3a='Max Likelihood: sigma=';
leg3b=' mean=';
leg3=strcat(leg3a,num2str(sigma),leg3b,num2str(x0));
  leg4=['1-sigma error for ' num2str(Nindep) ' indep. points'];
  [eN,dN,eA,dA]=binStatistics(bins,1,Nindep,4*pi*(180/pi)^2,0);

  scale=4*pi*(180/pi)^2/Nindep;
  mychi2=sum( (((N/scale)-eN).^2)./eN );
  Q=1-chi2cdf(mychi2,length(bins)-1);
  fprintf('Chi2: %f\tQ(%f,%d) = %f\n',mychi2,mychi2,length(bins)-1,Q);

%  plot(fine,zideal,'r',fine,z,'g',bins,eA+dA,'r-.',bins,eA-dA,'r-.');
  plot(fine,zideal,'k-',fine,z,'k-.');
%  legend('Data','Ideal gaussian (sigma=1 mean=0)',leg3,leg4);
fprintf('Legend:\n%s\n%s\n%s\n%s\n','Data','Ideal gaussian (sigma=1 mean=0)',leg3,leg4);
grid;
hold off;
axis([bmin,bmax,0,max(N)*1.3]);
xlabel('SNR');
ylabel('sky area (deg^2)');
% title('Histogram of SNR');
cmap=colormap;
cmap(1,:)=[1,1,1]*.5;
colormap(cmap);
