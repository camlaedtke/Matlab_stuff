% Thomas cyclically symmetric attractor

conds = struct('x_0',0.1,'y_0',-0.25,'z_0',0.15,...
               'linewidth',1.5,'pointsize',1.5,...
               'reltol',1e-5,'abstol',1e-6,'stats','off');
conds.tstart = 0;
conds.tend = 500;
conds.gstart = 0.0;
conds.gend = 0.35;
conds.numsamples = 3500;

gamma = linspace(conds.gstart,conds.gend, conds.numsamples+1);
rangegamma = length(gamma);
conds.gamma = conds.gstart;
eq = getthomaseqns(conds);
eqns = [eq.x, eq.y, eq.z];
sol = odenumeric3D(eqns, conds);
p = plot3(sol.y(:,1), sol.y(:,2), sol.y(:,3), ...
    'XDataSource','sol.y(:,1)',...
    'YDataSource', 'sol.y(:,2)',...
    'ZDataSource','sol.y(:,3)', 'linewidth', conds.linewidth);

% Figure Properties
fig = get(groot,'CurrentFigure');
fig.Position = [100 100 1400 1000];
% Textbox
dim = [0.73 0.8 0.2 0.1];
dim2 = [0.73 0.72 0.2 0.1];
str1 = ['Relative Tolerance: ',num2str(conds.reltol)];
str2 = ['Absolute Tolarance: ',num2str(conds.abstol)];
str3 = ['Simulated Time: ',num2str(conds.tend),' seconds'];
str4 = ['Gamma: ', num2str(conds.gamma)];
str = {str1,str2,str3};
text = annotation(fig,'textbox',dim,'String',str,'FitBoxToText','on');
text.FontSize = 18;
text2 = annotation(fig,'textbox',dim2,'String',str4,'FitBoxToText','on');
text2.FontSize = 18;
% Axis
title(['Thomas Cyclically Symmetric Attractor \gamma = ',...
       num2str(conds.gstart),' to \gamma = ',num2str(conds.gend)]);
ax1 = gca;
ax1.Position = [0.15 0.15 0.7 0.7];
ax1.Color = 'White';
% Axis - limits
ax1.XLim = [-10 10];
ax1.YLim = [-10 10];
ax1.ZLim = [-10 10];
% Axis - labels
ax1.XLabel.String = 'x';
ax1.YLabel.String = 'y';
ax1.ZLabel.String = 'z';
ax1.LabelFontSizeMultiplier = 2.0;
ax1.TitleFontSizeMultiplier = 2.0;

grid on
linkdata 'on'

F(rangegamma) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('ThomasEqnsVideo', 'MPEG-4');
v.FrameRate = 30;
open(v)
camorbit(-30,-20)
pb = CmdLineProgressBar('Progress... ');

for i = 1:rangegamma  
    conds.gamma = gamma(i);
    eq = getthomaseqns(conds);
    eqns = [eq.x, eq.y, eq.z];
    S = odenumeric3D(eqns, conds);
    sol.y = zeros(length(S.y(:,1)), length(S.y(:,2)), length(S.y(:,3)));
    sol.y(:,1) = S.y(:,1);
    sol.y(:,2) = S.y(:,2); 
    sol.y(:,3) = S.y(:,3);
    text2.String = ['Gamma: ', num2str(conds.gamma)];
    refreshdata
    drawnow   
    camorbit(0.05, 0)
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
    pb.print(i,rangegamma)
end
close(v)
movie(fig,F,1)

function equ = getthomaseqns(conds)
    syms x(t) y(t) z(t) b
    dx = diff(x, t);
    dy = diff(y, t);
    dz = diff(z, t);
    equ.x = dx == (sin(y) - conds.gamma*x);
    equ.y = dy == (sin(z) - conds.gamma*y);
    equ.z = dz == (sin(x) - conds.gamma*z);
end