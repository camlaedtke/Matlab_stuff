function [x, y] = odelogseparation(s1, s2, NumPoints)
    x1  = s1.Var1;
    y1 = s1.Var2;
    x2  = s2.Var1;
    y2 = s2.Var2;
    ini   = min(min(x1), min(x2)); 
    fin   = max(max(x1), max(x2));
    x = linspace(ini, fin, NumPoints); 
    % Gotta make sure no duplicate points in the input to 'interp1'
    [x1, index] = unique(x1);
    y1 = interp1(x1, y1(index), x1);
    [x2, index] = unique(x2);
    y2 = interp1(x2, y2(index), x2);
    y1X = interp1(x1, y1, x, 'linear','extrap');
    y2X = interp1(x2, y2, x, 'linear','extrap');
    yDiff = y1X - y2X;
    y = log10(abs(yDiff));
end