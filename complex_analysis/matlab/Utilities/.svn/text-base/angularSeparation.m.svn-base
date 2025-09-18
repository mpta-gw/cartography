function gamma=angularSeparation(ra1,decl1,ra2,decl2);
% calculates the angular separation across a great arc
% between the two locations
% result in degrees


rd=pi/180;
alpha=rd*(90-decl1);
beta =rd*transpose(90-decl2);
C    = pi/12*(ra1*ones(fliplr(size(ra2))) - ones(size(ra1))*transpose(ra2));
gamma=180/pi*acos( cos(alpha)*cos(beta) + (sin(alpha)*sin(beta)).*cos(C) );

