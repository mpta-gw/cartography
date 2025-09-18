% review_plotter.m
%
% This routine demonstrate the Spherical harmonics transform and plotting
% routines used for plotting the Spherical Harmonics Analysis output
% 
% Approximations of 2-dimentinal delta up to Lmax are created at the north
% pole and then rotated to sky locations theta, phi
% The code also demonstrates the normalization, both numerically and
% direct.
%
% see also: testSignalFrom Northpole.m
% This function proves that calGammaLM.m (which is used in the analysis)
% is consistent with Spherical harmonics transform and plotting routines
%
%           It calculates the glm using calGammaLM.m 
%           and calculates transpose(plm)*glm.data
%           for a few Lmax-approximatons of a delta function.
%           Then it compares the resulting spectrum against
%           the theoretical expectation for a true point source,
%           which is just the overlap reduction function integrand:
%                gamma.gamma0.*exp(i*2*pi*freq.*gamma.tau)
%           
%           Plotted is the normalized difference of the two spectra.
%              Note that the deviation at high frequencies is real
%              bacause there is a light difference between a true point
%              source and a Lmax-approximation of a point source.


addpath(genpath('../../'));

Lmax=20;
theta=[pi/4,pi/2,3/4*pi,pi/2];
phi=[0,pi/2,pi,3/2*pi];

pp=0;
close all;
for kk=1:length(theta)
    % Y00 = 1/sqrt(4*pi)
    Y00=getYlm(Lmax,0,0);

    % create a point at theta, phi
    pole=deltaAtNorthPole(Lmax);
    R=WignerRotation(Lmax,theta(kk));
    Z=zRotation(Lmax,phi(kk));
    p=Z*R*pole;
    pp=pp+p;
    
    % demonstrate normalization:
    % integral int dOmega p = sqrt(4*pi)*Y00'*p != 1
    integral = sqrt(4*pi)*Y00'*p;
    
    % check the integral numerically
    [map,ra,decl]=makemap(p,1);
    weight = cos(decl/180*pi)*(pi/180).^2;
    display(['theta :',num2str(theta(kk)*180/pi),'deg, phi :',num2str(phi(kk)*12/pi),'h']);
    display(['Direct integral (00 component): ',num2str(integral)]);
    display(['Numerical Integral:             ',num2str(sum(map.*weight))]);

    % demonstrate location
    plotPlmMap(p);
    title(['theta :',num2str(theta(kk)*180/pi),'deg, phi :',num2str(phi(kk)*12/pi),'h']);
end

plotPlmMap(pp);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run testSignalFromNorthpole
testSignalFromNorthpole;

figure(6);
title(['theta :',num2str(0),'deg, phi :',num2str(0),'h']);
figure(7);
title(['theta :',num2str(0),'deg, phi :',num2str(0),'h']);
figure(8);
title(['theta :',num2str(90),'deg, phi :',num2str(0),'h']);
figure(9);
title(['theta :',num2str(90),'deg, phi :',num2str(0),'h']);
figure(10);
title(['theta :',num2str(112.5),'deg, phi :',num2str(9),'h']);
figure(11);
title(['theta :',num2str(112.5),'deg, phi :',num2str(9),'h']);


for kk=1:11
    figure(kk);
    print('-dpdf');
end