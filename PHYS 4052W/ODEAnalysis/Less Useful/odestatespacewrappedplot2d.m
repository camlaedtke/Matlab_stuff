function S = odestatespacewrappedplot2d(sol, conds)
    x = linspace(conds.tstart,conds.tend,(conds.tend-...
        conds.tstart)*conds.numpoints +(conds.tend-...
        conds.tstart)/conds.numpoints);
    y = deval(sol,x)';
    S.x = wrapToPi(y(:,1));
    S.y = y(:,2);
    plot(S.x, S.y, '-r', 'LineWidth', conds.linewidth)
    title('State Space Wrapped from -pi to pi')
    xlabel('theta (radians)')
    ylabel('Dtheta (radians)')
end
