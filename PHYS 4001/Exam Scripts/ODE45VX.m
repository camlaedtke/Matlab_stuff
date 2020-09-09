function xdot = dummy(t,x)
 % define constants
 g = 9.8; % meters/s^2
 m  = 0.430; % kg
 c = 0.07; % meters
 r = c/(2*pi); % meters
 V = (4/3)*pi*r^2; %meters^3
 p = m/V; % kg/meters^3
 v_T = 25; %m/s
 b = (p*V*g)/(v_T^2); % kg
 
 x_1 = x(1);
 x_dot = x(2);
 
 
 xdot = zeros(size(x));
 xdot(1) = x_dot;
 xdot(2) = -(b*x_dot^2)/m;
 

end
