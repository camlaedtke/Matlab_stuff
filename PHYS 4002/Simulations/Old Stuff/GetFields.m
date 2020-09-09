function [V,Ex,Ey,Ez,xG,yG,zG] = GetFields(inputs)
    
    tic
    N = inputs.N; kC = inputs.kC; Q = inputs.Q;
    xC = inputs.xC; yC = inputs.yC; zC = inputs.zC;
    minRx = inputs.minRx; minRy = inputs.minRy; minRz = inputs.minRz;
    minX = inputs.minX; minY = inputs.minY; minZ = inputs.minZ;
    maxX = inputs.maxX; maxY = inputs.maxY; maxZ = inputs.maxZ;
    minR = inputs.minR; Esat = inputs.Esat; Vsat = inputs.Vsat;
    
    % [3D] region
    x = linspace(minX,maxX,N);
    y = linspace(minY,maxY,N);
    z = linspace(minZ,maxZ,N);
    [xG, yG, zG] = meshgrid(x,y,z);

    
   % fields
    V = zeros(N,N,N);
    Ex = zeros(N,N,N); 
    Ey = zeros(N,N,N);
    Ez = zeros(N,N,N);

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
   
   [Xszx, Xszy, Xszz] = size(xG);
   [Eszx, Eszy, Eszz] = size(Ex);
   [Vszx, Vszy, Vszz] = size(V);
   
   fprintf(['\n Size ... \n Grid: %.0f x %.0f x %.0f \n '...
                           '   E: %.0f x %.0f x %.0f \n '...
                           '   V: %.0f x %.0f x %.0f \n '],...
            Xszx,Xszy,Xszz,Eszx,Eszy,Eszz,Vszx,Vszy,Vszz)
   toc
end

