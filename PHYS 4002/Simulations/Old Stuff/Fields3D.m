clear all
close all
clc

N = 101; a = 0.01; nQ = 50; 
% Dimensions of region / saturation levels
minX = -1.5; minY = -1.5;  minZ = -0.5; 
maxX = 1.5;  maxY = 1.5;   maxZ = 0.5;
minR = 1e-6; minRx = 1e-6; minRy = 1e-6; minRz = 1e-6;

% Charge distribution
% [xC,yC,zC,nQ] = mySphere(nQ);

xC = linspace(-0.001,0.001,nQ);
yC = linspace(-1,1,nQ);
[xC,yC] = meshgrid(xC,yC);
points = [xC(:),yC(:)];
xC = points(:,1);
yC = points(:,2);
zC = zeros(length(xC),1);
nQ = length(xC);

Q0 = 5*(1/nQ)*1e-6; 
Q = Q0.*ones(1,nQ);
szq = length(Q);
fprintf('Number of Charges: %.0f \n',szq)

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

[X,Y,Z,index] = ReduceGrid(1,1,N-1,xG,yG,zG);
[V_r,E_xx,E_yy,E_zz] = ReduceFieldGrid(index,V,Ex,Ey,Ez);
[E_xx,E_yy,E_zz,E_norm] = NormalizeEField(E_xx,E_yy,E_zz);

V_slice = 34000; linewidth = 1.7; E_scale = 0.5;

close all
figure(1)   
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.4 0.2 0.7 0.7]); 

% Plot E field
hold on
h = quiver3(X,Y,Z,E_xx,E_yy,E_zz,'autoscalefactor',E_scale,'MaxHeadSize',0.5);
set(h,'color',[0 0 1],'linewidth',linewidth)

% Plot Charges
hold on
scatter3(-a+xC,-a+yC,-a+zC,10,'r','filled')

% Plot V surface
[faces,verts,colors] = isosurface(X,Y,Z,V_r,V_slice,E_norm);
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
      'FaceColor','interp','EdgeColor','interp',...
      'FaceAlpha',0.3,'EdgeAlpha',0.3)

xlabel('x  [m]'); ylabel('y  [m]'); zlabel('z  [m]')
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
  