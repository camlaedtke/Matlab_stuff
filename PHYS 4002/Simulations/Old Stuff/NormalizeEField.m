function [E_xx,E_yy,E_zz,E_norm] = NormalizeEField(E_xx,E_yy,E_zz)
    tic
    E_xx = normalize(E_xx,'norm',2);
    E_yy = normalize(E_yy,'norm',2);
    E_zz = normalize(E_zz,'norm',2);
    E_norm = sqrt(E_xx.^2+E_yy.^2+E_zz.^2);
    
    toc
end

