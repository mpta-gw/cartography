load('sumXFisher_S5flat')
load('bad_S5flat')
fisher_opt=fisher_opt-fisher_bad;
x_opt=x_opt-x_bad;

sigma00=subFisher(fisher_opt,0)^-0.5;           %uncertainty of monopole moment
%   5.3593e-49
p00=inv(subFisher(fisher_opt,0))*subX(x_opt,0); %monopole moment
%   -2.5456e-50

h0=0.72;                             %Hubbles constant
H0=100*(h0/3.08568E19);              %
omega=p00*(2*pi^2*100^3*sqrt(4*pi))/(3*H0^2);
%   -1.0905e-07
omega_sig=sigma00*(2*pi^2*100^3*sqrt(4*pi))/(3*H0^2);
%   2.2959e-06
%%% set all modified eigenvalues to the smallest non-modified value
%%% [invfisher_opt,pCovar]=reginv(fisher_opt);             %reg. fisher matrix
%%% set all modified eigenvalues to infinity
[invfisher_opt,pCovar]=reginv(fisher_opt,12,2/3);      %reg. fisher matrix
p_opt=invfisher_opt*x_opt;                             %clean map
res=1;
[sigma,sigmaPix,RA,DECL,dOmg,U]=getSigmaMap(pCovar,res);  %map of sigma
[map,ra,decl]=makemap(p_opt,res);                         %map of P
sn=[map,sigmaPix];

%Bayesian UL map----------------------------------------------------
confidence=0.90;
mu_min=0;
%mu_max=5e-44;
mu_max=67*6e-48; %prior = (sph_map/radio_map) * (S4 radiometer limits)
for ii=1:length(map)
  [mode, mu_UL, bayesfactor] = ...
    bayesianUL(map(ii), sigmaPix(ii), confidence, mu_min, mu_max);
  ul_map(ii)=mu_UL;
end
%-------------------------------------------------------------------

figure
plotHistogram(sn,360,181,400);
print('-depsc2','SNR_histo_S5flat.eps')
print('-djpeg','SNR_histo_S5flat.jpeg')

figure
[Cl,Cl_sig]=constructCl(p_opt,invfisher_opt);   %Cls measure angular structure
errorbar([0:length(Cl)-1],Cl,Cl_sig,'x')
set(gca,'FontSize',20); 
xlabel('l')
ylabel('Cl')
axis([0,20,-1e-96,1e-96])
print('-depsc2','Cl_S5flat.eps')
print('-djpeg','Cl_S5flat.jpeg')

save('S5flat.mat');
