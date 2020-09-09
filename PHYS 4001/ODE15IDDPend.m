function xdot = dummy(t, y0, yp0)
% Function arguments x = [theta(t), theta_dot(t)], t = [t initial, t final]

% Definte the phase state 
theta = y0;
theta_dot = yp0;

% Define a new vector:the derivative of our phase state
xdot = zeros(size(y0));
xdot(1) = theta_dot;
xdot(2) = pi*theta_dot*(-3/2) + ...
          pi^2*cos(t*pi*2)*(81/10)+ ...
          pi^2*(theta^3/6.0-theta)*9;