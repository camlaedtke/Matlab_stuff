syms w_x w_y w_z L
w = [w_x; w_y; w_z];
I = [1 0 1; 0 2 0; 1 0 1];
l = [1 0 0; 0 1 0; 0 0 1];

% getting eigenvalues
equ = det(I - L.*l) == 0
eigenvalues = solve(equ, L)

% getting principles axes
e_1 = (I - eigenvalues(1,:).*l)*w == [0; 0; 0]
e_2 = (I - eigenvalues(2,:).*l)*w == [0; 0; 0]
e_3 = (I - eigenvalues(3,:).*l)*w == [0; 0; 0]