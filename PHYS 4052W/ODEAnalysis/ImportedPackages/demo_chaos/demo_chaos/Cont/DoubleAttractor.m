function xs=DoubleAttractor(t,x,flag,rParam)

% System parameters

ms=730;
Ixx=2460;
w=6.9;
Fo=rParam;
ks1n=34000; 
ks2n=34000;
ks1=17500; 
ks2=17500;
coff1=784; 
coff2=784;
lf=1.011;
lr=1.803;
pha = rParam;%10

xg1=Fo*sin(w*t);
xg2=Fo*sin(w*t-pha);
xg1d=Fo*w*cos(w*t);
xg2d=Fo*w*cos(w*t-pha);

Xr=x(1)+lr*x(3);
Xf=x(1)-lf*x(3);
Xrd=x(2)+x(4)*lr;
Xfd=x(2)-x(4)*lf;

%state space vector
x1=x(1); x2=x(2); x3=x(3); x4=x(4);

xs=[x1 x2  x3  x4];

xs(1)=x2;

xs(2)=(-(coff2*Xrd+coff1*Xfd+(ks2*Xr+ks2n*Xr^3)+(ks1*Xf+ks1n*Xf^3))/ms)+((coff2*xg2d+coff1*xg1d+(ks2*xg2+ks2n*xg2^3)+(ks1*xg1+ks1n*xg1^3))/ms);

xs(3)=x4;

xs(4)=(-(coff2*Xrd*lr-coff1*Xfd*lf+(ks2*Xr+ks2n*Xr^3)*lr+(-ks1*Xf-ks1n*Xf^3)*lf)/Ixx)+((coff2*xg2d*lr-coff1*xg1d*lf+(ks2*xg2+ks2n*xg2^3)*lr+(-ks1*xg1-ks1n*xg1^3)*lf)/Ixx);

xs=xs';

end