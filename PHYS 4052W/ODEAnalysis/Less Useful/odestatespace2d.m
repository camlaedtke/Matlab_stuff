function odestatespaceplot2d(sol, conds)
    x = linspace(conds.tstart,conds.tend,(conds.tend-...
        conds.tstart)*conds.numpoints +(conds.tend-...
        conds.tstart)/conds.numpoints);
    y = deval(sol,x)';
    plot(y(:,1), y(:,2), '-r', 'LineWidth', conds.linewidth)
    title('State Space')
    xlabel('theta (radians)')
    ylabel('Dtheta (radians)')
end

