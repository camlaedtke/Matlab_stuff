% PHYS 4001 QUIZ 1\

Syms g k m q x1(t) x2(t);
assume(k>0);
assume(m>0);
assume(q>0);

D2x1 = diff(x1, t, 2);
D2x2 = diff(x1, t, 2);

ode1 = D2x1 == -2*(k/m)*x1 + (k/m)*x2 + 2*m*g*sin(q);
ode2 = D2x2 ==  (k/m)*x1 - (k/m)*x2 + m*g*sin(q);

odes = [ode1, ode2];
sol = dsolve(odes);
sol.x1
sol.x2


