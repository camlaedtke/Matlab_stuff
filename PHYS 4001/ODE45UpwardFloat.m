function ydot = dummy(t,y)
 g = 9.8; % m/s
 r = 1/3;
 V_T = 5;
 
 y_1 = y(1);
 y_dot = y(2);
 
 
 ydot = zeros(size(y));
 ydot(1) = y_dot;
 ydot(2) = g*(1-r)*(1 - y_dot^2/V_T^2);
 

end
