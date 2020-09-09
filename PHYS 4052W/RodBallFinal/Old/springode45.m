function S = springode45(N, w0, v_0, tspan)
    syms x(t) [N 1] 
    syms omega_0
    assume(omega_0 > 0)

    A = zeros(N,N);
    A(1,[1,2]) = [-1, 1];
    for i = 2:(N-1)
        A(i,[(i-1),i,(i+1)]) = [1, -2, 1];
    end 
    A(N,[N-1,N]) = [1, -1];
       
    DDx(t) = diff(x(t),t, 2);
    
    odes = DDx(t) == simplify(omega_0^2*A*x(t));
    [V] = odeToVectorField(odes);
%     if N == 384  
%         sym(S(617)) % x1
%         sym(S(618)) % dx1
%         sym(S(739))
%         sym(S(740)) % dx384        
%     end
    [V] = subs(V, omega_0, w0);
    M = matlabFunction(V, 'vars', {'t','Y'});
 
    Y_0 = zeros(1,2*N);
    
    % For large N, the locations of x, dx for the first and last
    % springs inside the function hangle are in wierd places. 
    % So we need to locate where the (dx1) element is in the vector
    % as well as the x(2N) and dx(2N-1) elements
    if N == 48 %
        Y_0(8) = v_0;
    elseif N == 96
        Y_0(104) = v_0;
    elseif N == 192
        Y_0(2*N-88) = v_0;
    elseif N == 384
        Y_0(618) = v_0;   
    else
        Y_0(4) = v_0;
    end
%     if N == 2
%         options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
%         'OutputFcn',@odephas2, 'OutputSel',[1 2]);
%     else
%         options = odeset('RelTol',1e-6,'AbsTol',1e-6,...
%         'OutputFcn',@odephas2, 'OutputSel',[2*N-1 2*N]);
%     end

%     options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    sol = ode45(M, tspan, Y_0);
    t_spr = sol.x;
    if N == 48
        du_spr = sol.y(5,:);
        v_spr = sol.y(6,:);
    elseif N == 96
        % sym(S(101))
        du_spr = sol.y(101,:);
        % sym(S(102))
        v_spr = sol.y(102,:);
    elseif N == 192
        du_spr = sol.y((2*N-1)-90,:);
        v_spr = sol.y(2*N-90,:);
    elseif N == 384
        du_spr = sol.y(739,:);
        v_spr = sol.y(740,:);
    elseif N == 2
        du_spr = sol.y(1,:);
        v_spr = sol.y(2,:);
    else
        du_spr = sol.y(2*N-1,:);
        v_spr = sol.y(2*N,:);
    end
    S = struct('t_spr',t_spr,'du_spr',du_spr,'v_spr',v_spr);
end