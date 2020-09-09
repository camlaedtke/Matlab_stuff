function S = odenumericsolve(equ, conds)
    [V] = odeToVectorField(equ);
    M = matlabFunction(V,'vars',{'t','Y'});
   
    Y_0 = [conds.y_0, conds.dy_0];
    options = odeset('RelTol', conds.reltol, 'AbsTol', conds.abstol, ...
              'Stats', conds.stats);
    tspan = [conds.tstart, conds.tend];
    S = ode45(M,tspan,Y_0, options);
end

