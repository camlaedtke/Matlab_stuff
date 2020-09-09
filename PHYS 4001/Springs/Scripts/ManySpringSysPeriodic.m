% normal modes by solving determinant = 0
% 8 objects & 8 springs with periodic boundary conditions 
% angular frequency = w, M/k = B

syms B w k m
assume(B>0)
assume(w>0)
% determinant is has non-trivial solutions only when EOM_matrix = 0
EOM_matrix = [-2+(m/k)*w^2 1 0 0 0 0 0 1;
              1 -2+(m/k)*w^2 1 0 0 0 0 0;
              0 1 -2+(m/k)*w^2 1 0 0 0 0;
              0 0 1 -2+(m/k)*w^2 1 0 0 0;
              0 0 0 1 -2+(m/k)*w^2 1 0 0;
              0 0 0 0 1 -2+(m/k)*w^2 1 0;
              0 0 0 0 0 1 -2+(m/k)*w^2 1;
              1 0 0 0 0 0 1 -2+(m/k)*w^2;];
          
% Make things cleaner              
M = subs(EOM_matrix, m/k, B);

equ = det(M) == 0;

% Solve the equation to get our normal modes, select only the unique ones
n_modes = unique(solve(equ, w))

syms A1 A2 A3 A4 A5 A6 A7 A8
assume(A1 ~=0);
V = [A1;A2;A3;A4;A5;A6;A7;A8];

M1 = subs(M, w, n_modes(1));
M2 = subs(M, w, n_modes(2));
M3 = subs(M, w, n_modes(3));
M4 = subs(M, w, n_modes(4));

eqns = M1*V == 0;
[A1,A2,A3,A4,A5,A6,A7,A8] = solve(eqns);
V1 = [A1,A2,A3,A4,A5,A6,A7,A8]';

eqns = M2*V == 0;
[A1,A2,A3,A4,A5,A6,A7,A8] = solve(eqns);
V2 = [A1,A2,A3,A4,A5,A6,A7,A8]';

eqns = M3*V == 0;
[A1,A2,A3,A4,A5,A6,A7,A8] = solve(eqns);
V3 = [A1,A2,A3,A4,A5,A6,A7,A8]';

eqns = M4*V == 0;
[A1,A2,A3,A4,A5,A6,A7,A8] = solve(eqns);
V4 = [A1,A2,A3,A4,A5,A6,A7,A8]';

e_vectors = [V1 V2 V3 V4]

% Test orthoganality
V1'*V2
V1'*V3
V1'*V4
V2'*V3
V2'*V4
V3'*V4