syms d1 d2 m_o m_c m_s u

D1 = 116*10^-12; % meters
D2 = 156*10^-12 % meters
M = 1.67*10^-27; % kg

% coordinates of our three objects
x = [-d1 0 d2];
y = [0 0 0];
z = [0 0 0];
% masses of our three objects
m = [m_o m_c m_s];


% Calculate moment of inertia tensor
Ixx =  ((m(1)*(y(1)^2 + z(1)^2)) + ...
        (m(2)*(y(2)^2 + z(2)^2)) + ...
        (m(3)*(y(3)^2 + z(3)^2)));
Ixy = -((m(1)*x(1)*y(1)) + ...
        (m(2)*x(2)*y(2)) + ...
        (m(3)*x(3)*y(3)));
Ixz = -((m(1)*x(1)*z(1)) + ...
        (m(2)*x(2)*z(2)) + ...
        (m(3)*x(3)*z(3)));
      
Iyx = -((m(1)*x(1)*y(1)) + ...
        (m(2)*x(2)*y(2)) + ...
        (m(3)*x(3)*y(3)));
Iyy =  ((m(1)*(x(1)^2 + z(1)^2)) + ...
        (m(2)*(x(2)^2 + z(2)^2)) + ...
        (m(3)*(x(3)^2 + z(3)^2)));
Iyz = -((m(1)*y(1)*z(1)) + ...
        (m(2)*y(2)*z(2)) + ...
        (m(3)*y(3)*z(3)));
    
Izx = -((m(1)*x(1)*z(1)) + ...
        (m(2)*x(2)*z(2)) + ...
        (m(3)*x(3)*z(3)));
Izy = -((m(1)*y(1)*z(1)) + ...
        (m(2)*y(2)*z(2)) + ...
        (m(3)*y(3)*z(3)));
Izz =  ((m(1)*(x(1)^2 + y(1)^2)) + ...
        (m(2)*(x(2)^2 + y(2)^2)) + ...
        (m(3)*(x(3)^2 + y(3)^2)));
    
I = [Ixx Ixy Ixz;
     Iyx Iyy Iyz;
     Izx Izy Izz];
 % I is in terms of mass*distance^2, which is correct
% Since the matrix is diagonal, it consists of the principle moments
I = simplify(I)
% substitute value of masses in terms of u
I = simplify(subs(I, [m_o, m_c, m_s], [16*u, 12*u, 32*u]))

% substitute value of u in terms of kg
% substitute values of d1 and d2 in terms of meters
I = subs(I, [u, d1, d2], [M, D1, D2]);
digits(3)
vpa(I)
% Answer makes sense. The molecule is confined to rotation 
% about y axis (Iyy) and about the z axis (Izz), both having the same
% value

% Since we are only interested in rotation about the z-axis, we
% can ignore the Iyy term

% With the COM as our origin
% Find distance from COM to each atom
% To do this, treat m_o and m_c as one combined mass,
% Use combined mass and m_s to calculate distance from origin the sulfer

% first we must get the location of our combined mass m_oc
syms m_oc f g  r
% m_oc is combined mass of oxigin and carbon
% Balancing the torques on each object...
equ = f*d1*m_o*g == (1-f)*d1*m_c*g;
f = simplify(solve(equ, f))

r_1 = (1-f)*d1;
r_2 = f*d1;

% r_1 is distance from carbon atom to m_oc
% r_2 is distance from oxygen atom to m_oc
r_1 = subs(r_1, [m_o, m_c], [16*u, 12*u])
r_2 = subs(r_2, [m_o, m_c], [16*u, 12*u]);

% Distance between sulfer atom and m_oc is d2 + r_1
% We will call that r


% calculate m_oc
% m_oc = m_o + m_c = 28u


r = d2 + r_1
r = subs(r, [d1, d2], [D1, D2]);
digits(3)
vpa(r) % = 222 pm, which is shorter than original sum of distances d1 + d2
% so that makes sense
r = d2 + r_1
% Now get the distance between the COM and m_oc and m_s
% Balancing the torque on both sides 

% distances of m_oc from COM is r_oc and distance m_s from COM is r_s
syms r_oc r_s
% Balancing torques
f = (m_s/(m_s + m_oc));

r_1 = (1-f)*r
r_2 = f*r

r_1 = subs(r_1, [m_oc, m_s], [28*u, 32*u])
r_2 = subs(r_2, [m_oc, m_s], [28*u, 32*u])

% distances from COM to m_oc and m_s are...
r_cm_oc = subs(r_1, [d1, d2], [D1, D2])
r_cm_s = subs(r_2, [d1, d2], [D1, D2])
digits(3)
vpa(r_1)
vpa(r_2)

% Finally, we can plug these distances into the moment of inertia tensor
% We just need Izz



% coordinates of our three objects
x = [0 -sym('r_cm_oc') sym('r_cm_s')];
y = [0 0 0];
z = [0 0 0];
% masses of our three objects
m = [0 m_oc m_s];


% Calculate moment of inertia tensor
Ixx =  ((m(1)*(y(1)^2 + z(1)^2)) + ...
        (m(2)*(y(2)^2 + z(2)^2)) + ...
        (m(3)*(y(3)^2 + z(3)^2)));
Ixy = -((m(1)*x(1)*y(1)) + ...
        (m(2)*x(2)*y(2)) + ...
        (m(3)*x(3)*y(3)));
Ixz = -((m(1)*x(1)*z(1)) + ...
        (m(2)*x(2)*z(2)) + ...
        (m(3)*x(3)*z(3)));
      
Iyx = -((m(1)*x(1)*y(1)) + ...
        (m(2)*x(2)*y(2)) + ...
        (m(3)*x(3)*y(3)));
Iyy =  ((m(1)*(x(1)^2 + z(1)^2)) + ...
        (m(2)*(x(2)^2 + z(2)^2)) + ...
        (m(3)*(x(3)^2 + z(3)^2)));
Iyz = -((m(1)*y(1)*z(1)) + ...
        (m(2)*y(2)*z(2)) + ...
        (m(3)*y(3)*z(3)));
    
Izx = -((m(1)*x(1)*z(1)) + ...
        (m(2)*x(2)*z(2)) + ...
        (m(3)*x(3)*z(3)));
Izy = -((m(1)*y(1)*z(1)) + ...
        (m(2)*y(2)*z(2)) + ...
        (m(3)*y(3)*z(3)));
Izz =  ((m(1)*(x(1)^2 + y(1)^2)) + ...
        (m(2)*(x(2)^2 + y(2)^2)) + ...
        (m(3)*(x(3)^2 + y(3)^2)));
    
I = [Ixx Ixy Ixz;
     Iyx Iyy Iyz;
     Izx Izy Izz];

 % Since the matrix is diagonal, it consists of the principle moments
I = vpa(simplify(I))
% I is in units of mass*distance^2, which is correct 
% substitute value of masses in terms of u
I = simplify(subs(I, [m_oc, m_s], [28*u, 32*u]))
% substitute value of distances in terms of meters
I = subs(I, [sym('r_cm_oc') sym('r_cm_s')], [r_cm_oc r_cm_s])
% substitute value of u in terms of kg
I = subs(I, u, M);
digits(3)
% Finally, our moment of inertia tensor is about the COM is
I_COM = vpa(I)

% Comparing the values of Izz we see there is less rotational inertia about
% the COM of the molecule than about the carbon atom, which makes sense

% Check both answers











