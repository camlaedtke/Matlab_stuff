% Number of grid point    [N = 1001]
   N = 21;
   a = 0.1;
   
   xC = [0,  0,  0, 0, 0];
   yC = [0,  0,  0, 0, 0];
    Q = [1, 0, 0, 0, 0] .* 1e-8;

   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);

   minX = -2;  maxX = 2;
   minY = -2;  maxY = 2;
   minZ = -2;  maxZ = 2;
   
   minR = 1e-6; minRx = 1e-6; minRy = 1e-6; minRz = 1e-6;

   Esat = 0.01* kC * max(abs(Q)) / a^2;
 
   V = zeros(N,N,N);
   Ex = zeros(N,N,N); Ey = zeros(N,N,N); Ez = zeros(N,N,N);
  
   x = linspace(minX,maxX,N);
   y = linspace(minY,maxY,N);
   z = linspace(minZ,maxZ,N);

   col1 = [1 0 0];
    
   [xG, yG, zG] = meshgrid(x,y,z);
    
for n = 1 : 5
   Rx = xG - xC(n);
   Ry = yG - yC(n);
   Rz = zG - zC(n);
   index = find(abs(Rx)+ abs(Ry)+ abs(Rz) == 0); 
   Rx(index) = minRx; 
   Ry(index) = minRy;
   Rz(index) = minRz;
   R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
   R(R==0) = minR;
   R3 = R.^3;
   Ex = Ex + kC .* Q(n) .* Rx ./ R3;
   Ey = Ey + kC .* Q(n) .* Ry ./ R3;
   Ez = Ez + kC .* Q(n) .* Rz ./ R3;
end

E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
if max(max(E)) >=  Esat; E(E >  Esat)  =  Esat; end;
if min(min(E)) <= -Esat; E(E < -Esat)  = -Esat; end;
   
if max(max(Ex)) >=  Esat; Ex(Ex >  Esat)  =  Esat; end;
if min(min(Ex)) <= -Esat; Ex(Ex < -Esat)  = -Esat; end;
   
if max(max(Ey)) >=  Esat; Ey(Ey >  Esat)  =  Esat; end;
if min(min(Ey)) <= -Esat; Ey(Ey < -Esat)  = -Esat; end;
   
if max(max(Ez)) >=  Esat; Ez(Ez >  Esat)  =  Esat; end;
if min(min(Ez)) <= -Esat; Ez(Ez < -Esat)  = -Esat; end;

close all
figure(1)    % 8888888888888888888888888888888888888888888888888888888888
set(gcf,'units','normalized','position',[0.4 0.1 0.6 0.8]); 
hold on
 index1 = 1 : 2 : N-1;
 index2 = index1;
 index3 = index1;
     
 X = xG(index1, index2, index3);
 Y = yG(index1, index2, index3);
 Z = zG(index1, index2, index3);
 
 % scaling of electric field lines: unit length
 E_x = Ex(index1, index2,index3)./(E(index1,index2,index3));
 E_y = Ey(index1, index2,index3)./(E(index1,index2,index3));
 E_z = Ez(index1, index2,index3)./(E(index1,index2,index3));
 
 % If these can be scaled in the right way, use it with colors
 
 E_xx = Ex(index1, index2, index3);
 E_yy = Ey(index1, index2, index3);
 E_zz = Ez(index1, index2, index3);
 
 E_xx = normalize(E_xx,'norm',1);
 E_yy = normalize(E_yy,'norm',1);
 E_zz = normalize(E_zz,'norm',1);
 
 R = sqrt(E_xx.^2+E_yy.^2+E_zz.^2);
 
 
 h = quiver3(X,Y,Z,E_xx,E_yy,E_zz,'autoscalefactor',2.8);
 
 set(h,'color',[0 0 1],'linewidth',0.8)
 
 box on
 view(-35,35)
 
 
 
 
 
 
 