close all
syms V r epsilon_0
assume(r>0)
q = 1;
V_0 = (4*pi*epsilon_0)/q;

x1 = 0;
y1 = r/sqrt(3);
z1 = 0;
R1 = sym(r*sqrt((x1+y1+z1)^2));

x2 = r/sqrt(3);
y2 = 2*r/sqrt(3);
z2 = 0;
R2 = sym(r*sqrt((x2+y2+z2)^2));

x3 = sym(r/sqrt(3));
y3 = sym(r/sqrt(3));
z3 = sym(r/sqrt(3));
R3 = sym(r*sqrt((x3+y3+z3)^2));

V = (1/(4*pi*epsilon_0))*(q/R1 + q/R2 - q/R3);
V_exact = V_0*V;

x = r/sqrt(3);
y = r/sqrt(3);
z = r/sqrt(3);

[phi, theta, r] = cart2sph(x,y,z);

V = V_0*(q/(4*pi*epsilon_0*r)+(q*cos(theta))/(2*pi*epsilon_0*r^2));
V_multi = simplify(V);


D = (V_multi - V_exact)/V_exact;

D = matlabFunction(D)

figure(1)
set(gcf,'units','normalized','position',[0.58 0.3 0.4 0.4],'MenuBar','none'); 
fplot(D,'LineWidth',3,'Color','b','DisplayName','Percent Diff')
title('$Percent Difference$','Interpreter','latex')
xlabel('$r$','Interpreter','latex')
ylabel('$Diff$','Interpreter','latex')
%xlim([0 2])
%ylim([0 10])
set(gca,'FontSize',20)
grid on


V_exact = matlabFunction(V_exact);
V_multi = matlabFunction(V_multi);


figure(2)
set(gcf,'units','normalized','position',[0.2 0.52 0.4 0.4],'MenuBar','none'); 
fplot(V_exact,'LineWidth',3,'Color','b','DisplayName','V_approx')
title('$V_{exact}$','Interpreter','latex')
xlabel('$r$','Interpreter','latex')
ylabel('$V_{0}$','Interpreter','latex')
xlim([0 2])
%ylim([0 10])
set(gca,'FontSize',20)
grid on

figure(3)
set(gcf,'units','normalized','position',[0.2 0.1 0.4 0.4],'MenuBar','none'); 
fplot(V_multi,'LineWidth',3,'Color','b','DisplayName','V_approx')
title('$V_{multi}$','Interpreter','latex')
xlabel('$r$','Interpreter','latex')
ylabel('$V_{0}$','Interpreter','latex')
xlim([0 2])
%ylim([0 10])
set(gca,'FontSize',20)
grid on







