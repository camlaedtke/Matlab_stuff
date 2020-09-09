%% PHYS 4002 HW5 PROBLEM 6
% Cameron Laedtke
% a) --------------------

% Potential V(x,y) for delta = 0.1a should start like a square wave
% but fall off on the edges.
% Potential for delta = a should have a higher maximum in the middle
% but smaller at the edges

close all
syms V_0 y n a delta
a = 1;
Vo = 1;
delta = 0.1*a;
V_0 = Vo*(delta/a)*tanh(y*(a-y)/delta^2);
V_0a = matlabFunction(V_0);

delta = a;
V_0 = Vo*(delta/a)*tanh(y*(a-y)/delta^2);
V_0b = matlabFunction(V_0);

figure(1) 
set(gcf,'units','normalized','position',[0.35 0.7 0.3 0.3]); 
p1 = fplot(V_0a,'LineWidth',2,'Color','b','DisplayName','delta = 0.1a');
hold on
p2 = fplot(V_0b,'LineWidth',2,'Color','r','DisplayName','delta = a');



% b) --------------------

% Keeping just the first five terms works well for delta = a because
% it looks "sine like", so lower order fourier series terms make it
% accurate, For delta = 0.1a, the potential looks more like a square wave,
% which needs higher order fourier terms to be well apporoximated (lots of
% error at the edges)

linesize = 1.5;
syms n a delta
a = 1;
Vo = 1;
x = 0;

delta = 0.1*a;
syms V y
V = 0;
fprintf("\nDelta = 0.1a \n")
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Va5 = matlabFunction(V);

delta = a;
syms V y
V = 0;
fprintf("\nDelta = a \n")
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Vb5 = matlabFunction(V);

p3 = fplot(Va5,'--','LineWidth',linesize,'Color','b',...
    'DisplayName','delta = 0.1a Est');
p4 = fplot(Vb5,'*','LineWidth',linesize,'Color','r',...
    'DisplayName','delta = a Est');
title('$V_{0}(y)$','Interpreter','latex')
xlabel('y/a','Interpreter','latex')
ylabel('$V_{0}$','Interpreter','latex')
xlim([0 1])
ylim([0 0.3])
set(gca,'FontSize',20)
legend([p1,p2,p3,p4],'NumColumns',2)
grid on
hold off

% Normalized error
err5a = @(y) abs((Va5(y) - V_0a(y)))./V_0a(y);
err5b = @(y) abs((Vb5(y) - V_0b(y)))./V_0b(y);

figure(3)
set(gcf,'units','normalized','position',[0.35 0.33 0.3 0.3]); 
fplot(err5a,'LineWidth',linesize,'Color','b')
hold on
fplot(err5b,'LineWidth',linesize,'Color','r')
title('Absolute Value of Normalized Difference','Interpreter','latex')
xlabel('y','Interpreter','latex')
ylabel('$|Norm  Diff|$','Interpreter','latex')
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',20)
grid on
hold off

% c)

% Potential decreases as x increases from zero

syms n a delta
a = 1;
Vo = 1;
x = 0.25;

delta = 0.1*a;
syms V y
V = 0;
fprintf("\nDelta = 0.1a \n")
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Va5x2 = matlabFunction(V);

delta = a;
syms V y
V = 0;
fprintf("\nDelta = a \n")
for n = 1:2:9
    Cn = getC(a,Vo,delta,n);
    V_n = Cn.*exp((-n.*pi.*x)./a).*sin((n.*pi.*y)./a);
    V = V+V_n; 
end
Vb5x2 = matlabFunction(V);

figure(4)
set(gcf,'units','normalized','position',[0.35 0.01 0.3 0.3]); 
p1 = fplot(Va5x2,'LineWidth',linesize,'Color','b',...
    'DisplayName','delta = 0.1a');
hold on
p2 = fplot(Vb5x2,'LineWidth',linesize,'Color','r',...
    'DisplayName','delta = a');

title('$V(a/4,y)$','Interpreter','latex')
xlabel('y/a','Interpreter','latex')
ylabel('$V_{0}$','Interpreter','latex')
xlim([0 1])
ylim([0 0.3])
set(gca,'FontSize',20)
legend([p1,p2])
grid on
hold off



function C = getC(a,Vo,delta,n)
    syms V_0(y)
    V_0(y) = Vo*(delta/a)*tanh(y*(a-y)/delta^2);
    F = V_0*sin((n*pi*y)/a);
    F = matlabFunction(F);
    y = (0:0.01:1)';
    C = (2/a)*trapz(y,F(y));
    fprintf('C(%.0f) = %.4f \n',n,C)
end