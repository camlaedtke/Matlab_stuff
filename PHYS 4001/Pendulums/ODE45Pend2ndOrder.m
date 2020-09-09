%dummy takes on the name of the .m script
function xdot = dummy(t,x)
g = 9.8; % m/s
R = 4; % m

% x is a 1x2 vector:
% x = [ theta    thetadot ]
theta=x(1);
thetadot=x(2);
% xdot is a 1x2 vector, same size as x
% xdot = [ thatadot   thetadoubledot]
xdot=zeros(size(x));
xdot(1)=thetadot;
% equation for 2nd order pendulum approx 
% thetadot(2)= - (g/R) (theta-(theta^3)/6) %g = 10 m/s^2; R = 0.5 m
xdot(2)= -(g/R)*(theta-(theta^3)/6);
end