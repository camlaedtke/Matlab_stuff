%// Create a quiver3 as we normally would (could also be 2D quiver)
clear all
close all
clc
tic


% INPUTS  ================================================================

% Number of grid point    [N = 1001]
   N = 201;
   
% Radius of circular charged conductor;   
   a = 0.05;
   
% X & Y components of position of charges 
    xCmin = -1; 
    xCmax = 1;
    yCmin = -1; 
    num_charges = 20;
    
    Q0 = (1/num_charges)*1e-6;
    Q = Q0.*ones(1,num_charges);
    NQ = length(Q);
    Q(1:NQ/2) = -Q0;
    
% Positions of charge distribution
    xC = linspace(xCmin,xCmax,num_charges);
    yC = xC.*zeros(1,num_charges);
    zC = xC.*zeros(1,num_charges);

% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
% Dimensions of region / saturation levels
%   [dimensions of region -2 to 2 / minR = 1e-6 / Esat = 1e6 / Vsat = 1e6]
   minX = -2;  maxX = 2;
   minY = -2;  maxY = 2;
   minZ = -2;  maxZ = 2;
   minR = 1e-6;
   minRx = 1e-6;
   minRy = 1e-6;
   minRz = 1e-6;
   Vsat = sqrt(num_charges/2) * kC * max(abs(Q)) / a;
   Esat = sqrt(num_charges/2) * kC * max(abs(Q)) / a^2;  
   fprintf('V saturation: %.2e \nE saturation: %.2e \n',Vsat,Esat)
   
% SETUP  =================================================================
   
  % fields
%     V = zeros(N,N);
%     Ex = zeros(N,N); 
%     Ey = zeros(N,N);
%     Ez = zeros(N,N);
    V = zeros(N,N,N);
    Ex = zeros(N,N,N); 
    Ey = zeros(N,N,N);
    Ez = zeros(N,N,N);
    
  % [2D] region
    x = linspace(minX,maxX,N);
    y = linspace(minY,maxY,N);
    z = linspace(minZ,maxZ,N);
    
  % color of charged object  +  red   /   - black
    col1 = [1 0 0];
    if Q(1) < 0; col1 = [0 0 0]; end;
        
  % grid positions
    [xG, yG, zG] = meshgrid(x,y,z);
    


% CALCULATION: POTENTIAL & ELECTRIC FIELD ================================

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

%%
close all
% GRAPHICS ===============================================================
     set(gcf,'units','normalized','position',[0.73 0.1 0.46 0.64]); 
     hold on
     index1 = 1 : 10 : N-1;
     %index1 = index1;
     index2 = index1;
     index3 = index1;
         
     X = xG(index1, index2, index3);
     Y = yG(index1, index2, index3);
     Z = zG(index1, index2, index3);
     
     % scaling of electric field lines: unit length
     U = Ex(index1, index2,index3)./(E(index1,index2,index3));
     V = Ey(index1, index2,index3)./(E(index1,index2,index3));
     W = Ez(index1, index2,index3)./(E(index1,index2,index3));


    q = quiver3(X, Y, Z, U, V, W);

    %// Compute the magnitude of the vectors
    mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), reshape(q.WData, numel(q.UData), [])).^2, 2));

    %// Get the current colormap
    currentColormap = colormap(gca);

    %// Now determine the color to make each arrow using a colormap
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));

    %// Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

    %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
    set(q.Head, ...
      'ColorBinding', 'interpolated', ...
      'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

    %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
    set(q.Tail, ...
      'ColorBinding', 'interpolated', ...
      'ColorData', reshape(cmap(1:2,:,:), [], 4).');






