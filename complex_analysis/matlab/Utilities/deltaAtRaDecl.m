function plm=deltaAtRaDecl(Lmax,RA,DECL)

% creates a delta function at RA (h) , DECL (deg)

phi=RA/12*pi;
theta=(90-DECL)/180*pi;
pole=deltaAtNorthPole(Lmax);
R=WignerRotation(Lmax,theta);
Z=zRotation(Lmax,phi);
plm=Z*R*pole;

