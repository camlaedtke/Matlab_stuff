%% Example 3.3

close all
a = 1;
V_0 = 1;
k = (4*V_0)/pi;
syms V x y
V = 0;
for n = 1:2:9
    V_0 = (1./n).*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_0; 
end
V = k.*V;
V = matlabFunction(V);

figure(1)
set(gcf,'units','normalized','position',[0.35 0.1 0.3 0.3]); 
fsurf(V,[0 1 0 1])
title('V(x,y)')
xlabel('x/a','Interpreter','latex')
ylabel('y/a','Interpreter','latex')
zlabel('$V/V_{0}$','Interpreter','latex')
set(gca,'FontSize',20)
view(30,20)

%exact solution
V_0 = 1;
a = 1;
V = @(x,y) ((2*V_0)/pi).*atan(sin((pi*y)/a)/sinh((pi*x)/a));

figure(2)
set(gcf,'units','normalized','position',[0.35 0.4 0.3 0.3]); 
fsurf(V,[0 1 0 1])
title('V(x,y)')
xlabel('x/a','Interpreter','latex')
ylabel('y/a','Interpreter','latex')
zlabel('$V/V_{0}$','Interpreter','latex')
zlim([0 1])
set(gca,'FontSize',20)
view(30,20)

%% Example 3.5
close all
a = 1;
b = 1;
k = 16/pi^2;
syms V y z
x = 1;
V = 0;
for n = 1:2:9
   for m = 1:2:9
   V_0 = (1./(n.*m)).*exp(-pi*sqrt((n./a).^2+(m./b).^2.*x)).*sin(((n.*pi.*y)./a)).*sin((m.*pi.*z)./b);
   V = V + V_0;
   end  
end

V = k.*V;
V = matlabFunction(V);
figure(1)
set(gcf,'units','normalized','position',[0.5 0.3 0.4 0.5]); 
fsurf(V,[0 1 0 1])
title('V(x,y)')
xlabel('y')
ylabel('z')
zlabel('V')
set(gca,'FontSize',20)
view(330,20)

%% PROBLEM 3.11
% Cameron Laedtke
close all
syms V x y 
z = 0.5;
a = 1;
k = 16/pi^2;
V = 0;
for n = 1:2:9
   for m = 1:2:9
   V_0 = (1/(n*m))*sin((n*pi.*x/a))*sin((m*pi*y/a))*((sinh(((pi*z)/a)*sqrt(n^2+m^2)))./sinh(pi.*sqrt(n^2+m^2)));
   V = V + V_0;
   end  
end

V = k.*V;
V = matlabFunction(V);

figure(1)
set(gcf,'units','normalized','position',[0.4 0.2 0.5 0.6]); 
fsurf(V,[0 1 0 1])
title('$V(x,y,a/2)$','Interpreter','latex')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('$V_{0}$','Interpreter','latex')
zlim([0 0.2])
set(gca,'FontSize',20)
view(330,20)
obtained = V(0.5,0.5);
expected = 1/6;
fprintf('\nExpected: %.4f Obtained: %.4f\n',expected,obtained)



