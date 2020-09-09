

conds = struct('tstart',0,'tend',50,...
               'y_0',-pi/2,'dy_0',0,...
               'linewidth',4,'pointsize',0.8,...
               'reltol',1e-7,'abstol',1e-8,'stats','off');
     
conds.gstart = 1.0;
conds.gend = 1.5;
conds.numsamples = 50;

conds.omega = 2*pi;
conds.omega_0 = 1.5*conds.omega;
conds.beta = conds.omega_0/4;

gamma = linspace(conds.gstart,conds.gend, conds.numsamples+1);
rangegamma = length(gamma);

conds.gamma = conds.gstart;
equ = getddpequation(conds);
sol = odenumericsolve(equ, conds);

p = plot3(sol.y(:,1),sol.y(:,2),sol.x(:,1), ...
    'XDataSource','sol.y(:,1)',...
    'YDataSource', 'sol.y(:,2)',...
    'ZDataSource','sol.x(:,1)', 'linewidth', conds.linewidth);

% Figure Properties
fig = get(groot,'CurrentFigure');
fig.Position = [100 100 1600 900];
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
title(['Damped Driven Pendulum \gamma = ',...
       num2str(conds.gstart),' to \gamma = ',num2str(conds.gend)]);
ax1 = gca;
ax1.Position = [0.1 0.18 0.85 0.75];
ax1.Color = 'White';
% Axis - limits
ax1.XLim = [-10 10];
ax1.YLim = [-5 5];
ax1.ZLim = [conds.tstart conds.tend];

% Axis - labels
ax1.XLabel.String = '\theta';
ax1.YLabel.String = '\partial \theta / \partial t';
ax1.ZLabel.String = 'time (seconds)';
ax1.LabelFontSizeMultiplier = 2.0;
ax1.TitleFontSizeMultiplier = 2.0;

grid on
linkdata 'on'

F(rangegamma) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('DDPVideo', 'MPEG-4');
v.FrameRate = 30;
open(v)
camorbit(0,-10)
pb = CmdLineProgressBar('Progress... ');

for i = 1:rangegamma  
    conds.gamma = gamma(i);
    equ = getddpequation(conds);
    S = odenumericsolve(equ, conds);
    sol.x = zeros(length(S.x(:,1)));
    sol.y = zeros(length(S.y(:,1)), length(S.y(:,2)));
    sol.x(:,1) = S.x(:,1);
    sol.y(:,1) = S.y(:,1); 
    sol.y(:,2) = S.y(:,2);
    text2.String = ['Gamma: ', num2str(conds.gamma)];
    refreshdata
    drawnow   
    camorbit(0.1, 0)
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
    pb.print(i,rangegamma)
end
close(v)
movie(fig,F,1)

function equ = getddpequation(conds)
    syms phi(t) 
    dphi = diff(phi, t);
    ddphi = diff(phi, t, 2);
    equ = ddphi + 2*conds.beta*dphi + conds.omega_0^2*sin(phi) == ...
    conds.gamma*conds.omega_0^2*cos(conds.omega*t);
end


function S = odenumericsolve(equ, conds)
    [V] = odeToVectorField(equ);
    M = matlabFunction(V,'vars',{'t','Y'});
    % Since it's a super stiff ode, we increase the tolerances
    Y_0 = [conds.y_0, conds.dy_0];
    options = odeset('RelTol', conds.reltol, 'AbsTol', conds.abstol, ...
              'Stats', conds.stats);
    tspan = [conds.tstart, conds.tend];
    S = ode45(M,tspan,Y_0, options);
    S.x = S.x';
    S.y = 0.1.*S.y';
end
