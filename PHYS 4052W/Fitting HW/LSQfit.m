function F = LSQfit(x,y,ey)
    sx = sum(x ./ (ey .^ 2) );
    sy = sum(y ./ (ey .^ 2) );
    sxx = sum((x .* x) ./ (ey .^ 2) );
    sxy = sum((x .* y) ./ (ey .^ 2) );
    s = sum(1 ./ (ey .^ 2) );
    delta=sxx*s-sx*sx;
    a=(sxx*sy-sx*sxy)/delta;
    ea=sqrt(sxx/delta);
    b=(s*sxy-sx*sy)/delta;
    eb=sqrt(s/delta);
    F =[ a, b ; ea, eb ];
end