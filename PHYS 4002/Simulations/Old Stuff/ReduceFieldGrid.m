function [V_r,E_xx,E_yy,E_zz] = ReduceFieldGrid(index,V,Ex,Ey,Ez)
    tic
    V_r = V(index, index, index);
    E_xx = Ex(index, index, index);
    E_yy = Ey(index, index, index);
    E_zz = Ez(index, index, index);
    
    [Eszx, Eszy, Eszz] = size(E_xx);
    [Vszx, Vszy, Vszz] = size(V_r);
    
     fprintf(['\n Size \n E: %.0f x %.0f x %.0f \n '...
                         'V: %.0f x %.0f x %.0f \n '],...
              Eszx,Eszy,Eszz,Vszx,Vszy,Vszz);
    
    toc
    
    
end

