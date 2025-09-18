function [eN,dN,eA,dA]=binStatistics(bins,sigm,N,A,plotFlag)
defval('sigm',1);
defval('N',1);
defval('A',1);
defval('plotFlag',2);

nB=length(bins);
binEdges=zeros(1,nB+1);
binEdges(2:nB)=1/2*(bins(1:nB-1)+bins(2:nB));
binEdges(1)=2*bins(1)-binEdges(1);
binEdges(nB+1)=2*bins(nB)-binEdges(nB);

bEerf=erf(binEdges./sigm./sqrt(2))./2;
p=bEerf(2:nB+1)-bEerf(1:nB);

eN=N.*p;
dN=sqrt(N.*p.*(1-p));
eA=A./N.*eN;
dA=A./N.*dN;
if or(plotFlag==1,plotFlag==3)
  figure;
  plot(bins,eN,bins,eN+dN,bins,eN-dN);
  grid;
  legend('<N>','<N>+dN','<N>-dN');
  xlabel('bins');
  ylabel('N');
end
if or(plotFlag==2,plotFlag==3)
  figure;
  plot(bins,eA,bins,eA+dA,bins,eA-dA);
  grid;
  legend('<A>','<A>+dA','<A>-dA');
  xlabel('bins');
  ylabel('A');
end;
return;
