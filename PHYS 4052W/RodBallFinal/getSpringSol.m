function [t,u,v] = getSpringSol(N,kball,tspan)
    load('tek0039ImportantParams.mat','m_ball','m_rod','w0')
    load('tek0039PendulumRemoved.mat','v_cm')
    beta = 3.334; % from fit to envelope peaks
    wball = sqrt(kball/m_ball);
    % v0ball = N*v_cm/2; % Changes particle velocity and CM velocity
    v0ball = (v_cm*(m_ball+m_rod))./(2*m_ball);
    y0ball = 0;
    w1 = sqrt(w0^2 - beta^2);
    w0_n = (w1*N)/pi;

    Y_0 = zeros(1,2*(N+1));
    Y_0(2*N+1) = y0ball;
    Y_0(2*N+2) = v0ball;

    options = odeset('stats','on','RelTol',1e-9,'AbsTol',1e-10); 
    S = ode45(@(t,y)springdde(t,y,w0_n,wball,N,beta), tspan, Y_0, options);

    t = S.x';
    t = t*10^3;
    u_left = S.y(1,:)';
    u = S.y((2*N-1),:)';
    u=u*10^3;
    v = S.y(2*N,:)';
    u_ball = S.y(2*N+1,:)';
    [~, ~, ~] = contactdata(u_ball, u_left, t, wball);
end

