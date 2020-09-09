function xdot = dummy(t, x)
% Function arguments x = [theta(t), theta_dot(t)],
% t = [t initial, t final]

% Definte the phase state 
theta = x(1);
theta_dot = x(2);

% Define a new vector:the derivative of our phase state
xdot = zeros(size(x));
xdot(1) = theta_dot;
xdot(2) = pi*theta_dot*(-3/2) + ...
          pi^2*cos(t*pi*2)*(81/10)+ ...
          pi^2*(theta^3/6.0-theta)*9;