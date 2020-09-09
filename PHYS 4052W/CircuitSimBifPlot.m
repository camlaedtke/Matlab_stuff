% CHAOTIC CIRCUIT
clear all
close all

% time = 100 - 400 should be good
time = 80; % seconds

V_0 = 0.25; % V
R_0 = 157000; % Ohms
R = 47000; 
C = 1*10^-6; % Farads
R_1 = 15000; 
R_2 = 88900; 
R_nu = 50000:100:110000; % Chaos between ~65k and ~110k
% paper has 2000 step resolution, which is 60000:35:120000

gamma = R_2/R_1; % controls type of chaotic behavior
omega = R/R_0; % controls amplitude of signal
beta = R./R_nu; % controls chaos
rangebeta = length(beta);

freq_0= 1/(2*pi*R*C); % Characteristic frequency
tau = 2*pi*R*C;
t_f = round(time/(R*C));% dimensionless time

% Define Initial Conditions and Plot Parameters
Y_0 = [V_0,0,0]; % first val is V0 (0.25V from paper)

tstart = 0;
tend = t_f;
bifurcationstart = 10; % (seconds)
pointsize = 5;
  
% For the plotting, I start with a blank scatterplot, and specify x and y
% data that will be constantly updated
% This is necessary because you can't just store everything
% in one huge array and plot it all at once, 

clf
% Matlab likes it when you pre-allocate memory by specifying size of arrays
timerange = tend - bifurcationstart;
D.x = zeros(timerange, rangebeta);
D.y = zeros(timerange, rangebeta);
t = linspace(0,timerange,tend);
V = linspace(0,timerange,tend);

f1 = figure(1);
f1.Units = 'normalized';
set(gcf,'position',[0.31 0.1 0.35 0.85],'MenuBar','none');
scatter(D.x(:), D.y(:), pointsize,...
       'XDataSource','D.x(:)','YDataSource', 'D.y(:)');
xlim([R_nu(1) R_nu(length(R_nu))])
ylim([-0.5 3.5])
title('Circuit Bifurcation Plot','interpreter','latex')
xlabel('Resistance of Variable Resistor (Ohms)','interpreter','latex')
ylabel('Local Maxima (V)','interpreter','latex')
set(gca,'FontSize',22,'Position',[0.12 0.11 0.84 0.84])
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
tspan = [0,tend];
grid on
linkdata 'on'


f2 = figure(2);
set(gcf,'units','normalized','position',[0.7 0.1 0.35 0.85],'MenuBar','none');
scatter(t,V,5,'red','XDataSource','t','YDataSource', 'V')
title('Voltage vs Time','interpreter','latex')
xlabel('Time','interpreter','latex')
ylabel('Voltage','interpreter','latex')
ax2 = gca;
ylim([-0.5 3.5])
set(gca,'FontSize',22,'Position',[0.1 0.11 0.84 0.84])
grid on
linkdata 'on' % tells the scatter plot to update its data

for i = 1:rangebeta  
    beta_0 = beta(i); % new resistor value
    % update function hangle with new resistor value
    M = @(t,Y)[Y(2);Y(3);...
          (-omega - beta_0.*Y(3)-gamma.*Y(1) - Y(2))*(Y(1)<0) + ...
          (-omega - beta_0.*Y(3)-Y(2))*(Y(1)>0)];
    S = ode45(M,tspan,Y_0,options);  % get solution
    t = S.x'; 
    y = S.y'; 
    V = y(:,1);
    [pks,locs] = findpeaks(V,t); % get peaks
    if i == 1
        ax2.XLim = [locs(50) locs(70)];
        ax2.XTickMode = 'manual';
        ax2.XTick = [];
    end
    refreshdata(f2)
    drawnow
    D.x(1:length(pks(150:end)),i) = R_nu(i);
    D.y(1:length(pks(150:end)),i) = pks(150:end); 
    refreshdata(f1)
    drawnow 
    pause(0.01);
end


