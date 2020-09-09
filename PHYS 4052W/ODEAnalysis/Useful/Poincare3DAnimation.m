%% Poincare Section Animation

%%
% Create 3-D scatter plot animation that shows how the Poincare section of 
% a second order differential equation changes with respect to a parameter

%%
% *Define Parameters of Differential Equation*

params = struct('omega',(2*pi));
params.omega_0 = 1.5*params.omega;
params.beta = params.omega_0/8;
%%
% *Define Initial Conditions and Plot Parameters*

conds = struct('tstart',0,'tend',1300,...
               'y_0',-pi/2,'dy_0',0,...
               'bifurcationstart',500,'poincarestart',100,...
               'linewidth',0.1,'pointsize',0.8,...
               'reltol',1e-5,'abstol',1e-5,'stats','off');
conds.gstart = 0.975;
conds.gend = 1.025;
conds.numpoints = conds.tend - conds.poincarestart;

numsamples = 5000;
gamma = linspace(conds.gstart,conds.gend, numsamples+1);
rangegamma = length(gamma);

clf
D.y = zeros(conds.numpoints, conds.numpoints, conds.numpoints);
%%
% The parameter |gamma| is a vector of gamma values from
% |conds.gstart| to |conds.gend|
%%
% The struct field |D.y| contains coordinates of our state space
% (and scatter plot) for a single value of |gamma| over an interval of 
% time defined by |conds.numpoints|. To pre-allocate memory,
% its size is defined. 

%%
% *Define Figure and Axes Properties*

p = scatter3(D.y(:,1), D.y(:,2),D.y(:,3), conds.pointsize,...
    'XDataSource','D.y(:,1)',...
    'YDataSource', 'D.y(:,2)',...
    'ZDataSource','D.y(:,3)');

fig = get(groot,'CurrentFigure');
%%
% |fig.Position| coordinates 1 and 2 are the distance from bottom left of screen.
% Coordinates 3 and 4 are the width and height of the figure (in pixels)

fig.Position = [100 100 1800 1000];
ax1 = gca;
%%
% |ax1.Position| defines the amount of space between the edges of the axes
% to the edges of the figure window 

ax1.Position = [0.15 0.15 0.7 0.7];
ax1.Color = 'White';
ax1.LabelFontSizeMultiplier = 2.0;
ax1.TitleFontSizeMultiplier = 2.0;
ax1.XLim = [-3.2 3.2];
ax1.YLim = [-10 30];
ax1.ZLim = [conds.gstart conds.gend];
ax1.XAxisLocation = 'bottom';
ax1.XLabel.String = '\theta';
ax1.YLabel.String = '\partial \theta / \partial t';
ax1.ZLabel.String = '\gamma';
title(['Poincare Section from \gamma = ',num2str(conds.gstart),...
    ' to \gamma = ',num2str(conds.gend)]);
linkdata 'on'

%%
% *Animation*
% Each frame is stored in frame object |F|. Pre-allocate space in |F|.
% |VideoWriter| is the parent object for all animation objects

F(rangegamma) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('PoincareSectionVideo2', 'MPEG-4');
v.FrameRate = 60;
open(v)
camorbit(-30,-20)
pb = CmdLineProgressBar('Progress... ');
%%
% *Iterate Over Gamma Values*
%%
% Perform the following operations at each iteration:
%%
%
% # Update parameter value |params.gamma|
% # Redefine the starting time as |conds.tstart| = 0
% # Get the differential equation with updated parameters
% # Solve the differential equation numerically with updated parameters
% # Redefine the starting time |conds.tstart = conds.poincarestart|
% # Get the Poincare section for the updated solution
% # Replace the old plotted data with the new data
% # Link the new data to the plot
% # Generate and save frame
%

for i = 1:rangegamma  
    params.gamma = gamma(i);
    conds.tstart = 0;
    equ = odegetddpequation(params); 
    sol = odenumericsolve(equ, conds); 
    conds.tstart = conds.poincarestart; 
    S = odepoincaredata(sol, conds);  
    D.y(:,1) = S.x;
    D.y(:,2) = S.y; 
    D.y(:,3) = params.gamma*ones(size(S.y));
    refreshdata
    drawnow   
    camorbit(0.01,0.002)
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
    pb.print(i,rangegamma)
end
close(v)
%movie(fig,F,1)

function S = odenumericsolve(equ, conds)
    [V] = odeToVectorField(equ);
    M = matlabFunction(V,'vars',{'t','Y'});
   
    Y_0 = [conds.y_0, conds.dy_0];
    options = odeset('RelTol', conds.reltol, 'AbsTol', conds.abstol, ...
              'Stats', conds.stats);
    tspan = [conds.tstart, conds.tend];
    S = ode45(M,tspan,Y_0, options);
end

function S = odepoincaredata(sol, conds)
    % Poincare plot, plots a state space point each cycle
    t_cycle = linspace(conds.tstart,conds.tend-1,(conds.tend-conds.tstart));
    y_cycle = deval(sol,t_cycle)';
    x_cycle = wrapToPi(y_cycle(:,1));
    S.x = x_cycle;
    S.y = y_cycle(:, 2);
end

function equ = odegetddpequation(params)
    syms phi(t) 
    dphi = diff(phi, t);
    ddphi = diff(phi, t, 2);
    equ = ddphi + 2*params.beta*dphi + params.omega_0^2*sin(phi) == ...
    params.gamma*params.omega_0^2*cos(params.omega*t);
end


