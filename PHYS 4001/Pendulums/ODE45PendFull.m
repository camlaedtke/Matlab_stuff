% This function defines the differential equations being solved
function ydot=dummy(t,y)
% y is a 2x1 column vector:
% y = [ theta    ]
%     [ thetadot ]
theta=y(1);
thetadot=y(2);
% doty is a 2x1 column vector:
% y = [ thetadot   ]
%     [ thetadoubledot ]
ydot=zeros(2, 1);
ydot(1)=thetadot;
% full equation for pendulum, take g/length = 20
g = 9.8; % m/s
R = 4; % m
ydot(2)=-(g/R)*sin(theta);
end