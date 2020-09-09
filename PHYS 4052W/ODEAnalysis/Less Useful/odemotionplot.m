function S = odemotionplot(sol, conds)
    x = linspace(conds.tstart,conds.tend,(conds.tend-...
        conds.tstart)*conds.numpoints +(conds.tend-...
        conds.tstart)/conds.numpoints);
    y = deval(sol,x)';
    S.x = x';
    S.y = y(:,1);
    plot(S.x, S.y, '-r', 'LineWidth', conds.linewidth)
    title('Motion Plot')
    xlabel('time (seconds)')
    ylabel('theta (radians)')
end

