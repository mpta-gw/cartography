function RZ=centerRaDecl(Lmax,RA,DECL)

% Rotates by RA (h) , DECL (deg)

phi=(12-RA)/12*pi;
theta=(0-DECL)/180*pi;

R=WignerRotation(Lmax,theta);
Z=zRotation(Lmax,phi);
RZ=R*Z;

