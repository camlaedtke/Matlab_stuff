%% Numerical integration
close all
linesize = 1.5;
syms n a delta
a = 1;
Vo = 1;
% x = 0;


delta = 0.1*a;
syms V x y
V = 0;
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Va5 = matlabFunction(V);

delta = a;
syms V x y
V = 0;
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Vb5 = matlabFunction(V);

figure(1)
set(gcf,'units','normalized','position',[0.35 0.43 0.3 0.3]); 
fsurf(Va5)
title('V(0,y)','Interpreter','latex')
xlabel('x/a','Interpreter','latex')
ylabel('y/a','Interpreter','latex')
zlabel('$V/V_{0}$','Interpreter','latex')
xlim([0 1])
ylim([0 1])
zlim([0 0.3])
set(gca,'FontSize',20)
grid on
view(30,20)

figure(2)
set(gcf,'units','normalized','position',[0.35 0.11 0.3 0.3]); 
fsurf(Vb5)
title('V(0,y)','Interpreter','latex')
xlabel('x/a','Interpreter','latex')
ylabel('y/a','Interpreter','latex')
zlabel('$V/V_{0}$','Interpreter','latex')
xlim([0 1])
ylim([0 1])
zlim([0 0.3])
set(gca,'FontSize',20)
grid on
view(30,20)



function C = getC(a,Vo,delta,n)
    syms V_0(y)
    V_0(y) = Vo*(delta/a)*tanh(y*(a-y)/delta^2);
    F = V_0*sin((n*pi*y)/a);
    F = matlabFunction(F);
    y = (0:0.01:1)';
    C = (2/a)*trapz(y,F(y));
    %fprintf('C(%.0f) = %.7f \n',n,C)
end