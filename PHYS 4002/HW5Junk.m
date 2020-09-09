close all
syms x

P0 = legendreP(0,cos(x));
P1 = legendreP(1,cos(x));
P2 = legendreP(2,cos(x));
P3 = legendreP(3,cos(x));
P4 = legendreP(4,cos(x));
P5 = legendreP(5,cos(x));
P6 = legendreP(6,cos(x));
P7 = legendreP(7,cos(x));
P8 = legendreP(8,cos(x));
P9 = legendreP(9,cos(x));
P10 = legendreP(10,cos(x));
P11 = legendreP(11,cos(x));
P12 = legendreP(12,cos(x));
P13 = legendreP(13,cos(x));

figure(1)
set(gcf,'units','normalized','position',[0.3 0.5 0.65 0.5]); 
fplot(P0)
hold on
fplot(P1)
fplot(P2)
fplot(P3)
fplot(P4)
fplot(P5)
fplot(P6)
fplot(P7)
fplot(P8)
fplot(P9)
fplot(P10)
fplot(P11)
fplot(P12)
fplot(P13)
title('Legendre Polynomial in cos(x)')
xlim([0 2*pi])
ylim([-1 1])
grid on
hold off

P0 = legendreP(0,cos(x));
P1 = P0 + legendreP(1,cos(x));
P2 = P1 + legendreP(2,cos(x));
P3 = P2 + legendreP(3,cos(x));
P4 = P3 + legendreP(4,cos(x));
P5 = P4 + legendreP(5,cos(x));
P6 = P5 + legendreP(6,cos(x));
P7 = P6 + legendreP(7,cos(x));
P8 = P7 + legendreP(8,cos(x));
P9 = P8 + legendreP(9,cos(x));
P10 = P9 + legendreP(10,cos(x));
P11 = P10 + legendreP(11,cos(x));
P12 = P11 + legendreP(12,cos(x));
P13 = P12 + legendreP(13,cos(x));

figure(2)
set(gcf,'units','normalized','position',[0 0.25 0.65 0.5]); 
fplot(P0)
hold on
fplot(P1)
fplot(P2)
fplot(P3)
fplot(P4)
fplot(P5)
fplot(P6)
fplot(P7)
fplot(P8)
fplot(P9)
fplot(P10)
fplot(P11)
fplot(P12)
fplot(P13)
title('Legendre Polynomial in cos(x) Summed')
xlim([0 2*pi])
ylim([0 14])
grid on
hold off

figure(3)
set(gcf,'units','normalized','position',[0.3 0.01 0.65 0.5]); 
fplot(int(legendreP(0,cos(x))))
hold on
fplot(int(legendreP(1,cos(x))))
fplot(int(legendreP(2,cos(x))))
fplot(int(legendreP(3,cos(x))))
fplot(int(legendreP(4,cos(x))))
fplot(int(legendreP(5,cos(x))))
fplot(int(legendreP(6,cos(x))))
fplot(int(legendreP(7,cos(x))))
fplot(int(legendreP(8,cos(x))))
fplot(int(legendreP(9,cos(x))))
fplot(int(legendreP(10,cos(x))))
fplot(int(legendreP(11,cos(x))))
fplot(int(legendreP(12,cos(x))))
fplot(int(legendreP(13,cos(x))))
title('Legendre Polynomial in cos(x) Integrals')
xlim([0 2*pi])
ylim([-pi pi])
grid on
hold off

