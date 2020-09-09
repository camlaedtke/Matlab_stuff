function S = odenumeric3D(eqns, conds)
    [V] = odeToVectorField(eqns);
    M = matlabFunction(V,'vars',{'t','Y'});
    % Since it's a super stiff ode, we increase the tolerances
    Y_0 = [conds.x_0, conds.y_0, conds.z_0];
    options = odeset('RelTol', conds.reltol, 'AbsTol', conds.abstol, ...
              'Stats', conds.stats);
    %tspan = conds.tstart:conds.sep:conds.tend
    S = ode45(M,conds.tspan,Y_0, options);
    S.x = S.x';
    S.y = S.y';
end
