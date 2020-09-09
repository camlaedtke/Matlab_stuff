function xdot = dummy(t, x)
% Function arguments x = [theta(t), theta_dot(t)], 
                    %t = [t initial, t final]

% Set initial conditions
g = 9.8; % meters/second^2
L = 2; % meters
mu = 0.3;

% Definte the phase state 
theta = x(1);
theta_dot = x(2);

% Define a new vector:the derivative of our phase state
% Define xdot = [theta_dot, theta_double_dot]
xdot = zeros(size(x));
xdot(1) = theta_dot;
xdot(2) = -mu*theta_dot - (g/L)*sin(theta);
% The vector [xdot(1), xdot(2)] is now the input to the next iteration
