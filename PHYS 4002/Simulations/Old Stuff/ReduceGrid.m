function [X, Y, Z, index] = ReduceGrid(start,step,stop,xG,yG,zG)
    tic
    index = start : step : stop; 

    X = xG(index, index, index);
    Y = yG(index, index, index);
    Z = zG(index, index, index);
    
    [Xszx, Xszy, Xszz] = size(X);
    fprintf('\n Size ... \n Grid: %.0f x %.0f x %.0f \n ',Xszx,Xszy,Xszz);
    toc
end

