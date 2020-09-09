function S = odestatespaceplot2d(sol, conds)
    x = linspace(conds.tstart,conds.tend,(conds.tend-...
        conds.tstart)*conds.numpoints +(conds.tend-...
        conds.tstart)/conds.numpoints);
    y = deval(sol,x)';
    S.x = y(:,1);
    S.y = y(:,2);
    plot(S.x, S.y, '-r', 'LineWidth', conds.linewidth)
    title('State Space')
    xlabel('theta (radians)')
    ylabel('Dtheta (radians)')
  
end

