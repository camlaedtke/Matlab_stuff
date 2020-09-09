%% Define parameters, initial conditions
clear all 
close all

params = struct('omega',(2*pi));
params.gamma = 1.5;
params.omega_0 = 1.5*params.omega;
params.beta = params.omega_0/8;
conds = struct('tstart',10,'bifurcationstart',500,'tend',300,...
    'y_0',-pi/2,'dy_0',0,'linewidth',0.1,'pointsize',3,...
    'reltol',1e-5,'abstol',1e-5,'stats','off');
conds.numpoints = conds.tend - conds.bifurcationstart;
conds.gstart = 01.070;
conds.gend = 1.087;
conds.samplefreq = 10000; 
conds.poincarestart = 10;

%% Poincare Section
conds.numpoints = conds.tend - conds.poincarestart;
% range of gamma values to sample
gamma = linspace(conds.gstart,conds.gend,...
         (conds.gend-conds.gstart)*conds.samplefreq+1);
rangegamma = length(gamma);

% State space [theta, Dtheta] represented as D.y(:,1), D.y(:,2)
% For each value of gamma, there will be a state space matrix
% that describes how [theta, Dtheta] evolves over time

D.y = zeros(conds.numpoints, conds.numpoints);
p = scatter(D.y(:,1), D.y(:,2), conds.pointsize,...
    'XDataSource','D.y(:,1)','YDataSource', 'D.y(:,2)');


fig = get(groot,'CurrentFigure');
fig.Position = [200 200 1000 1000];
ax1 = gca;
ax1.Position = [0.1 0.1 0.87 0.82];
ax1.Color = 'White';
ax1.LabelFontSizeMultiplier = 1.5;
ax1.TitleFontSizeMultiplier = 1.75;
ax1.XLim = [-3.5 3.5];
ax1.YLim = [-20 40];
ax1.XAxisLocation = 'bottom';
ax1.XLabel.String = 'Theta';
ax1.YLabel.String = 'DTheta';
title('Poincare Diagram')
linkdata 'on'

F(rangegamma) = struct('cdata',[],'colormap',[]);
v = VideoWriter('PoincareSectionVideoTest', 'MPEG-4');
open(v)

% TODO: add an extra culumn for gamma value associated with 
% each state space point
pb = CmdLineProgressBar('Progress ...');
for i = 1:rangegamma
    params.gamma = gamma(i);
    conds.tstart = 0;
    equ = odegetddpequation(params);
    sol = odenumericsolve(equ, conds);
    conds.tstart = conds.poincarestart;
    S = odepoincaredata(sol, conds);
    
    D.y(:,1) = S.x;
    D.y(:,2) = S.y;  
    
    refreshdata
    drawnow
    
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
    pb.print(i,rangegamma)
end
close(v)
movie(fig,F,2)

