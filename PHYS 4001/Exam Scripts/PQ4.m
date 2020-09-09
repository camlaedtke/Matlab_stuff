% Using ODE45 for position
theta = pi/4 % radians
v_0 = 50 % m/s

x_o = 0;
dx_o = v_0*cos(theta);
y_o = 0;
dy_o = v_0*sin(theta);


tspan = [0, 5]
X_0 = [x_o, dx_o]
Y_0 = [y_o, dy_o]

[t,x] = ode45(@ODE45VX, tspan, X_0)
[t,y] = ode45(@ODE45VY, tspan, Y_0)

clf

plot(t,y(:,1), '*r')
xlabel('Time (s)')
ylabel('Range (m)')

plot(t,x(:,1), '*r')
xlabel('Time (s)')
ylabel('Range (m)')
% we see the time 