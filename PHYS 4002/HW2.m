%% HOMEWORK 2 PROBLEM 7
close all
clear all

figure(1) 
set(gcf,'units','normalized','position',[0.55 0.01 0.4 0.5]); 

syms x y L f(x)
y = 12/16;
f = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f,'DisplayName','12L/16')
title('V(x,y)')
xlabel('x distance from rod')
ylabel('Potential V (V_0)')
xlim([0 5])
ylim([0 6])
set(gca,'FontSize',20)
hold on

y = 11/16;
f2 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f2,'DisplayName','11L/16')

y = 10/16;
f3 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f3,'DisplayName','10L/16')

y = 9/16;
f4 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f4,'DisplayName','9L/16')

y = 8/16;
f5 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f5,'DisplayName','8L/16')

y = 7/16;
f6 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f6,'DisplayName','7L/16')

y = 6/16;
f7 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f7,'DisplayName','6L/16')

y = 5/16;
f8 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f8,'DisplayName','5L/16')

y = 4/16;
f9 = log((sqrt(x.^2+y.^2)+y)./(sqrt(x.^2+y.^2)-y));
fplot(f9,'DisplayName','4L/16')
legend()
grid on
hold off

figure(2) 
syms x y L f(x,y)
set(gcf,'units','normalized','position',[0.55 0.5 0.4 0.5]); 
y = L;
f = log((sqrt(x^2+y^2)+y)/(sqrt(x^2+y^2)-y));
fsurf(f)
xlim([0 1])
ylim([0 5])
zlim([0 6])
zlabel('Potential V (V_0)')
title('V(x,y) for y = 0 to y = L')
xlabel('L')
ylabel('x distance from rod')
set(gca,'FontSize',20)
view(140,20)

