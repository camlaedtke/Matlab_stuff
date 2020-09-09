function equ = getddpequation(params)
    syms phi(t) 
    dphi = diff(phi, t);
    ddphi = diff(phi, t, 2);
    equ = ddphi + 2*params.beta*dphi + params.omega_0^2*sin(phi) == ...
    params.gamma*params.omega_0^2*cos(params.omega*t);
end

