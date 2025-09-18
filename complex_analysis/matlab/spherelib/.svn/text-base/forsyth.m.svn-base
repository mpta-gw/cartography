function varargout=forsyth(Te,lbd,f,rc,drho,T);
% [G2b,k,l,Zb,Zf]=FORSYTH(Te,lbd,f,rc,drho,T);
%
% Calculates the predicted coherence-square between Bouguer anomalies and
% topography over a range of wavenumbers for a lithosphere loaded at the
% surface and at one other depth which is also the depth of
% compensation. It doesn't matter where this depth is; it is implicitly
% included in the loading factor.
%
% Works for a vector of Te OR a vector of f (but not both).
%
% INPUT:
%
% Te      elastic thickness of plate [m]
% lbd     the wavelength [m]
% f       ratio of bottom-to-top applied loads
%         0 for surface loading only (produces unity)
%         1 for equal load
%         infinity for moho loading only (produces unity)
% rc      density contrast across topography interface [kg/m^3]
%         rho_crust for loading of continents
% drho    density contrast at compensation interface
%         (usually = rho_mantle - rho_crust) [kg/m^3]
% T       depth to the density contrast/compensation interface, in m
%
% OUTPUT:
%
% G2b       Bouguer coherence-square (between 0 and 1)
% k         wavenumber [rad/m]
% l         wavelength [m]
% Zb        Bouguer addmittance, in mgal/m
%
% See also XTRAXIS1D, MCKENZIE
%
% EXAMPLE:
%
% forsyth('demo')
%
% Last modified by fjsimons-at-alum.mit.edu, October 21st, 2003.

if ~isstr(Te) & ~strcmp(Te,'demo')
  defval('Te',[20 80]*1e3);
  defval('lbd',linspace(10e3,2000e3,100))
  defval('f',1);
  defval('rc',2670);
  defval('drho',630);
  defval('T',40e3)
  defval('E',1.4e11);
  v=0.25;
  g=9.81;
  G=fralmanac('GravCst');
  disp(sprintf('E= %5.3g; v= %5.3f',E,v))

  k=2*pi./lbd;

  % Turcotte and Schubert (3-72)
  D=(E*Te.^3)/(12*(1-v^2)); % Flexural Rigidity [Pa*m^3 or N*m]

  % Create grid on which to calculate coherence-square
  if length(Te)>=1 & length(f)==1
    [LL,DD]=meshgrid(lbd,D);
    FF=f;
  elseif length(f)>=1 & length(Te)==1
    [LL,FF]=meshgrid(lbd,f);
    DD=D;
  else
    error('Only one vector allowed')
  end
  KK4=(2*pi./LL).^4;
  KK=(2*pi./LL);

  % Forsyth Eqs. (3) and (6)
  xai=1+DD.*KK4/drho/g;
  phi=1+DD.*KK4/rc/g;
  beta=rc./xai/drho;

  % See Forsyth Eq. (25)
  Ctop=(xai.*drho^2+FF.^2.*rc^2.*phi).^2;
  Cbot1=xai.^2.*drho^2+FF.^2.*rc^2;
  Cbot2=drho^2+FF.^2.*rc^2.*phi.^2;
  G2b=Ctop./Cbot1./Cbot2;

  % Bouguer admittance
  Zb=-2*pi*G*rc*exp(-KK*T).*(phi.*FF.^2.*beta.^2+1./xai)./...
		   (FF.^2.*beta.^2+1);
  % Free air admittance
  Zf=Zb+2*pi*G*rc;
  
  % Convert to mgal/m
  Zb=Zb/1e-5;
  Zf=Zf/1e-5;

  l=lbd;
  
  varnames={'G2b','k','l','Zb','Zf'};
  for ind=1:nargout
    varargout{ind}=eval(varnames{ind});
  end
else
  % Illustrates the functions FORSYTH, MCKENZIE, and TRANSL
  clf
  [G2bF,k,l,ZbF,ZfF]=forsyth([20 80]*1e3,[],1,2670,630);
  subplot(211)
  pF=semilogx(k*1000,G2bF,'Color','g','LineW',2);
  [k12,l12]=transl(1,[20 80],2670,630);
  [l,Zb20,G2b20,Zf20]=mckenzie([0 2670 2670+630],[0 40],[1 1],20);
  [l,Zb80,G2b80,Zf80]=mckenzie([0 2670 2670+630],[0 40],[1 1],80);
  hold on
  pM20=semilogx(2*pi./l*1000,G2b20,'b');
  pM80=semilogx(2*pi./l*1000,G2b80,'r');
  set([pM20 pM80],'MarkerS',4); grid on; openup(gca,6); longticks(gca);
  xlim([3e-3 1e-1])
  xl(1)=xlabel('Wavenumber (rad/km)');
  yl(1)=ylabel('Bouguer Coherence \gamma^2');
  x1=xtraxis1d(gca); longticks(x1)
  hold on
  pl(1)=plot(l12(1),0.5,'x');
  pl(2)=plot(l12(2),0.5,'o');
  xl(2)=xlabel('Wavelength (km)');
  subplot(212)
  pZ=semilogx(k*1000,ZfF*1000); hold on
  pZ20=semilogx(2*pi./l*1000,Zf20*1000,'b');
  pZ80=semilogx(2*pi./l*1000,Zf80*1000,'r');
  grid on; openup(gca,6); longticks(gca);
  xlim([3e-3 1e-1])
  xl(3)=xlabel('Wavenumber (rad/km)');
  yl(2)=ylabel('Free-air Admittance Z_f (mgal/km)');
  x2=xtraxis1d(gca); longticks(x2)
  figdisp
end
