close all
% create a 2D grid
th = linspace(0,pi,400);    % inclination
phi = linspace(0,2*pi,400); % azimuth
[th,phi] = meshgrid(th,phi);

% compute spherical harmonic of degree 3 and order 1
Y1 = harmonicY(12,6,th,phi,'type','real');
Y2 = harmonicY(14,4,th,phi,'type','real');
Y3 = harmonicY(12,7,th,phi,'type','real');
% plot the magnitude
r1 = abs(Y1);
r2 = abs(Y2);
r3 = abs(Y3);
r = r1+r2+r3;


[x,y,z] = sph2cart(phi,pi/2-th,r);
figure(1)
set(gcf,'units','normalized','position',[0.37 0 0.63 0.9]); 
surf(x,y,z,r);
title('$Spherical Harmonics$','Interpreter','latex')
set(gca,'FontSize',20)
