function S = odestatespacedata2d(sol, conds)
%     x = linspace(conds.tstart,conds.tend,(conds.tend-...
%         conds.tstart)*conds.numpoints +(conds.tend-...
%         conds.tstart)/conds.numpoints)
    x = linspace(conds.tstart,conds.tend-1, conds.numpoints);
    y = deval(sol,x)';
    S.x = y(:,1);
    S.y = y(:,2);
end

