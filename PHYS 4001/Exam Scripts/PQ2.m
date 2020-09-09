syms theta(t) phi(t) psi(t) k m l omega

v1 = l*diff(theta);
v2 = l*diff(phi);
v3 = l*diff(psi);

h1 = l*(1-cos(theta));
h2 = l*(1-cos(phi));
h3 = l*(1-cos(psi));

U_k1 = (1/2)*k*l^2*(sin(theta)-sin(phi))^2;
U_k2 = (1/2)*k*l^2*(sin(phi)-sin(psi))^2;

T = (1/2)*m*(v1^2 + v2^2 + v3^2);
U = m*g*h1 + m*g*h2 + m*g*h3 + U_k1 + U_k2;

L = T-U;
% Euler Lagrance equation w.r.t each generalized coordinate and velocity
EOM_theta = functionalDerivative(L,theta) == 0;
EOM_phi = functionalDerivative(L,phi) == 0;
EOM_psi = functionalDerivative(L,psi) == 0;

% use small angle approximation
old = [sin(theta), sin(phi), sin(psi), cos(theta), cos(phi), cos(psi)];
new = [theta, phi, psi, (1-(theta^2)/2), (1-(phi^2)/2), (1-(psi^2)/2)];

EOM_theta = subs(EOM_theta, old, new);
EOM_phi = subs(EOM_phi, old, new);
EOM_psi = subs(EOM_psi, old, new);
EOM_theta = simplify(EOM_theta);
EOM_phi = simplify(EOM_phi);
EOM_psi = simplify(EOM_psi);

M = [EOM_theta;
     EOM_phi;
     EOM_psi];

syms A B C
assume(l>0)
assume(m>0)

theta_o = A*exp(i*omega*t);
dtheta = diff(theta_o);
ddtheta = diff(dtheta);

phi_o = B*exp(i*omega*t);
dphi = diff(phi_o);
ddphi = diff(dphi);

psi_o = C*exp(i*omega*t);
dpsi = diff(psi_o);
ddpsi = diff(dpsi);
old = [theta, diff(theta, t), diff(theta, t, 2), ...
       phi, diff(phi, t), diff(phi, t, 2),...
       psi, diff(psi, t), diff(psi, t, 2)];
new = [theta_o, dtheta, ddtheta,...
       phi_o, dphi, ddphi,...
       psi_o, dpsi, ddpsi];
   
M = subs(M, old, new);
old = [exp(omega*i*t),exp(2*omega*i*t)];
new = [1, 1];

M = subs(M, old, new);
M = simplify(M)
