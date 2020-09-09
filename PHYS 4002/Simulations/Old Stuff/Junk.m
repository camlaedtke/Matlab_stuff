clear all
close all
clc

N = 101;
a = 0.05;
nQ = 50; % num charges
Q0 = 5*(1/nQ)*1e-6; 
Q = Q0.*ones(1,nQ);

% Dimensions of region / saturation levels
minX = -2;  maxX = 2;
minY = -2;  maxY = 2;
minZ = -2;  maxZ = 2;
minR = 1e-6; minRx = 1e-6; minRy = 1e-6; minRz = 1e-6;

% Positions of charge distribution
xCmin = -1; xCmax = 1; % charge start/end
xC = linspace(xCmin,xCmax,nQ);
yC = sin(2.*xC).*ones(1,nQ);
zC = xC.*zeros(1,nQ);

% Constants / Saturation levels
eps0 = 8.854e-12; kC = 1/(4*pi*eps0);
Esat = sqrt(nQ/2) * kC * max(abs(Q)) / a^2; 
Vsat = sqrt(nQ/2) * kC * max(abs(Q)) / a;

% color of charged object  (+  red  /  - black)
col1 = [1 0 0]; col2 = [0 0 0];

inputs = struct('N',N,'a',a,'eps0',eps0,'kC',kC,'Esat',Esat,'Vsat',Vsat,...
'Q',Q,'nQ',nQ,'xC',xC,'yC',yC,'zC',zC,'col1',col1,'col2',col2,...
'minX',minX,'minY',minY,'minZ',minZ,'maxX',maxX,'maxY',maxY,'maxZ',maxZ,...
'minR',minR,'minRx',minRx,'minRy',minRy,'minRz',minRz);


[V,Ex,Ey,Ez,xG,yG,zG] = GetFields(inputs);

%% Graph

[X,Y,Z,index] = ReduceGrid(1, 3 ,N-1,xG,yG,zG);
[V_r,E_xx,E_yy,E_zz] = ReduceFieldGrid(index,V,Ex,Ey,Ez);
[E_xx,E_yy,E_zz,E_norm] = NormalizeEField(E_xx,E_yy,E_zz);

V_slice = 35000; linewidth = 0.6; E_scale = 1.0;

close all
figure(1)   
set(gcf,'units','normalized','position',[0.1 0.6 0.5 0.6]); 
hold on
h = quiver3(X,Y,Z,E_xx,E_yy,E_zz,'autoscalefactor',E_scale,'MaxHeadSize',0.5);
set(h,'color',[0 0 1],'linewidth',linewidth)

hold on
for n = 1:nQ
    pos = [-a+xC(n), -a+yC(n), 2*a, 2*a];
    h = rectangle('Position',pos,'Curvature',[1,1]);
    set(h,'FaceColor',col1,'EdgeColor',col1);
end
   
[faces,verts,colors] = isosurface(X,Y,Z,V_r,V_slice,E_norm);
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
'FaceColor','interp','EdgeColor','interp',...
'FaceAlpha',0.3,'EdgeAlpha',0.3)

xlabel('x  [m]');
ylabel('y  [m]');
zlabel('z  [m]')
title('|L2 NORM| direction of scaled E at grid points',...
 'fontweight','normal');
set(gca,'xLim',[minX,maxX]); 
set(gca,'yLim',[minY,maxY]);
set(gca,'zLim',[minZ,maxZ]);
axis([minX maxX minY maxY minZ maxZ]);

view(35,20)
axis equal
grid on
box on
lighting gouraud
  