function ydot = dummy(t,y)
 % define constants
 g = 9.8; % meters/s^2
 m  = 0.430; % kg
 c = 0.07; % meters
 r = c/(2*pi); % meters
 V = (4/3)*pi*r^2; %meters^3
 p = m/V; % kg/meters^3
 v_T = 25; %m/s
 b = (p*V*g)/(v_T^2); % kg
 
 y_1 = y(1);
 y_dot = y(2);
 
 
 ydot = zeros(size(y));
 ydot(1) = y_dot;
 ydot(2) = -(b*y_dot^2)/m - g;
 

end