clear all
close all
clc
tic

N = 201;
a = 0.05;

linecharge = true;

% constants
eps0 = 8.854e-12;
kC = 1/(4*pi*eps0);

% Dimensions of region / saturation levels
minX = -2; minY = -2; minZ = -2;
maxX =  2; maxY =  2; maxZ =  2;
minR = 1e-6;
minRx = 1e-6; minRy = 1e-6; minRz = 1e-6;

% fields
V = zeros(N,N,N);
Ex = zeros(N,N,N); 
Ey = zeros(N,N,N);
Ez = zeros(N,N,N);

% [3D] region
x = linspace(minX,maxX,N);
y = linspace(minY,maxY,N);
z = linspace(minZ,maxZ,N);
[xG, yG, zG] = meshgrid(x,y,z);

% color of charged object  +  red   /   - black
col1 = [1 0 0]; col2 = [0 0 0];

if linecharge
   xCmin = -1; xCmax = 1; nQ = 50;
   Q0 = 5*(1/nQ)*1e-6;
   Q = Q0.*ones(1,nQ);
   Q(1:nQ/2) = -Q0;
   % Positions of charge distribution
   xC = linspace(xCmin,xCmax,nQ);
   yC = sin(2.*xC).*ones(1,nQ);
   zC = xC.*zeros(1,nQ);
   Vsat = sqrt(nQ/2) * kC * max(abs(Q)) / a;
   Esat = sqrt(nQ/2) * kC * max(abs(Q)) / a^2; 
   for n = 1 : length(Q)
        Rx = xG - xC(n);
        Ry = yG - yC(n);
        Rz = zG - zC(n);

        index = find(abs(Rx) + abs(Ry) + abs(Rz) == 0);
        Rx(index) = minRx;  
        Ry(index) = minRy;
        Rz(index) = minRz;

        R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
        R(R==0) = minR;
        V = V + kC .* Q(n) ./ (R);

        R3 = R.^3;
        Ex = Ex + kC .* Q(n) .* Rx ./ R3;
        Ey = Ey + kC .* Q(n) .* Ry ./ R3;
        Ez = Ez + kC .* Q(n) .* Rz ./ R3;
   end
else
   Q = [20, 0, 0, 0, 0] .* 1e-6;
   xC = [0,  0,  0, 0, 0];
   yC = [0,  0,  0, 0, 0];
   zC = [0,  0,  0, 0, 0];
   Vsat = kC * max(abs(Q)) / a;
   Esat = kC * max(abs(Q)) / a^2;
   col1 = [1 0 0];
   for n = 1 : 5
        Rx = xG - xC(n);
        Ry = yG - yC(n);
        Rz = zG - zC(n);

        index = find(abs(Rx)+ abs(Ry)+abs(Rz) == 0); 
        Rx(index) = minRx;  
        Ry(index) = minRy;
        Rz(index) = minRz;

        R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
        R(R==0) = minR;
        V = V + kC .* Q(n) ./ (R);

        R3 = R.^3;
        Ex = Ex + kC .* Q(n) .* Rx ./ R3;
        Ey = Ey + kC .* Q(n) .* Ry ./ R3;
        Ez = Ez + kC .* Q(n) .* Rz ./ R3;
   end
end

if max(max(V)) >=  Vsat; V(V > Vsat)  = Vsat; end;
if min(min(V)) <= -Vsat; V(V < -Vsat) = -Vsat; end;

E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
if max(max(E)) >=  Esat; E(E >  Esat)  =  Esat; end;
if min(min(E)) <= -Esat; E(E < -Esat)  = -Esat; end;

if max(max(Ex)) >=  Esat; Ex(Ex >  Esat)  =  Esat; end;
if min(min(Ex)) <= -Esat; Ex(Ex < -Esat)  = -Esat; end;

if max(max(Ey)) >=  Esat; Ey(Ey >  Esat)  =  Esat; end;
if min(min(Ey)) <= -Esat; Ey(Ey < -Esat)  = -Esat; end;

if max(max(Ez)) >=  Esat; Ez(Ez >  Esat)  =  Esat; end;
if min(min(Ez)) <= -Esat; Ez(Ez < -Esat)  = -Esat; end;

toc
  %%
close all

V_slice = 30000; %
linewidth = 0.6;
vector_scale = 1.6;

figure(1)   
 set(gcf,'units','normalized','position',[0.1 0.6 0.3 0.4]); 
 hold on
 index1 = 11 : 15 : N-11; index2 = index1; index3 = index1;

 X = xG(index1, index2, index3);
 Y = yG(index1, index2, index3);
 Z = zG(index1, index2, index3);
 
 V_r = V(index1, index2, index3);
 E_xx = Ex(index1, index2, index3);
 E_yy = Ey(index1, index2, index3);
 E_zz = Ez(index1, index2, index3);

 E_normalize = 2;
 E_xx = normalize(E_xx,'norm',E_normalize);
 E_yy = normalize(E_yy,'norm',E_normalize);
 E_zz = normalize(E_zz,'norm',E_normalize);

 R_norm = sqrt(E_xx.^2+E_yy.^2+E_zz.^2);

 h = quiver3(X,Y,Z,E_xx,E_yy,E_zz,...
     'autoscalefactor',vector_scale,'MaxHeadSize',0.5);
 set(h,'color',[0 0 1],'linewidth',linewidth)

 hold on
 % charges   
 if linecharge
    col = col2;
    for n = 1:nQ/2
        pos = [-a+xC(n), -a+yC(n), 2*a, 2*a];
        h = rectangle('Position',pos,'Curvature',[1,1]);
        set(h,'FaceColor',col,'EdgeColor',col);
    end
    col = col1;
    for n = 1+nQ/2:nQ
        pos = [-a+xC(n), -a+yC(n), 2*a, 2*a];
        h = rectangle('Position',pos,'Curvature',[1,1]);
        set(h,'FaceColor',col,'EdgeColor',col);
    end
 else
    pos = [-a+xC(n), -a+yC(n), 2*a, 2*a];
    h = rectangle('Position',pos,'Curvature',[1,1]);
    set(h,'FaceColor',col1,'EdgeColor',col1);
 end

 xlabel('x  [m]');
 ylabel('y  [m]');
 zlabel('z  [m]')
 title('|L2 NORM| direction of scaled E at grid points',...
     'fontweight','normal');

 set(gca,'xLim',[minX,maxX]); 
 set(gca,'yLim',[minY,maxY]);
 set(gca,'zLim',[minZ,maxZ]);
 axis([minX maxX minY maxY minZ maxZ]);
 axis equal
 grid on
 box on

 [Rx, Ry, Rz] = size(R_norm);
 [Eszx, Eszy, Eszz] = size(E_xx);
 [Xszx, Xszy, Xszz] = size(X);
 [Vszx, Vszy, Vszz] = size(V_r);

 fprintf(['\n Size \n R: %.0f x %.0f x %.0f \n '...
                   'E: %.0f x %.0f x %.0f \n '...
                   'X: %.0f x %.0f x %.0f \n '...
                   'V: %.0f x %.0f x %.0f \n '],...
                   Rx,Ry,Rz,Eszx,Eszy,Eszz,Xszx,Xszy,Xszz,...
                   Vszx,Vszy,Vszz)

 [faces,verts,colors] = isosurface(X,Y,Z,V_r,V_slice,R_norm);
 patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
'FaceColor','interp','EdgeColor','interp',...
'FaceAlpha',0.3,'EdgeAlpha',0.3)
 view(35,20)

 lighting gouraud
  

 %%
     toc