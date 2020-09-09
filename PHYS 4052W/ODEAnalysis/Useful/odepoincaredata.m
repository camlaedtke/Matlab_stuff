function S = odepoincaredata(sol, conds)
    % Poincare plot, plots a state space point each cycle
    t_cycle = linspace(conds.tstart,conds.tend-1,(conds.tend-conds.tstart));
    y_cycle = deval(sol,t_cycle)';
    x_cycle = wrapToPi(y_cycle(:,1));
    S.x = x_cycle;
    S.y = y_cycle(:, 2);
end

