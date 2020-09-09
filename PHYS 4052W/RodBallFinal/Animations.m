clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')
% load('100MassFirst3ms.mat')
% Get vector of 100 colors ranging from blue to red
r = (0:0.01:1)';
g = zeros(101,1);
b = (1:-0.01:0)';
rgb = [r,g,b];

figure_size = [0.42 0.05 0.56 0.8];

%%

close all

load('tek0039Sim100M10msecTc119usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~, fb, idx_contact] = loadSol(S,N,wball);

t_range = [0.00 tau]; 
xlim_shift = -tau; 
font_size = 15;

F(48) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRod2DTc119msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

% Plot of displacement vs time for all 100 masses simultaniously
close all
figure(1)
subplot(2,1,1)
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
plot(t_sim(idx_contact), 0.4*fb*10^-11,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot(t_sim,S.y((i*2-1),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize',font_size)
ax1 = gca;
ax1.XLim = t_range+xlim_shift;
xlim(t_range+xlim_shift) 

subplot(2,1,2)
% Plot of acceleration vs time for all 100 masses simultaniously
hold on
for i=1:100
    plot(t_sim(1:end-1),diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ax2 = gca;
ax2.XLim = t_range+xlim_shift;

% Animation parameters
dt = tau/48;

% Animate first period or so
for j = 1:48
    set(ax1, 'XLim', ax1.XLim+dt)
    set(ax2, 'XLim', ax2.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)




%% Animate
close all

load('tek0039Sim100M10msecTc119usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

% Define video object 

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc119msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.1])
zlim([-2500 2500])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

% Animation parameters
dt = tau/480;

% Animate first period or so
for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)
%

close all
load('tek0039Sim100M10msecTc110usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc110msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.12])
zlim([-3000 3000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)
%

close all
load('tek0039Sim100M10msecTc100usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

% Define video object 

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc100msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.12])
zlim([-4000 4000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)
%

close all
load('tek0039Sim100M10msecTc090usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

% Define video object 

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc090msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.14])
zlim([-5000 4000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)
%
close all
load('tek0039Sim100M10msecTc080usLowTol.mat')

[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

% Define video object 

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc080msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.16])
zlim([-6000 6000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)

%
close all
load('tek0039Sim100M10msecTc070usLowTol.mat')
[t_sim, ~, ~, ~, ~, ~, ~] = loadSol(S,N);

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

% Define video object 

F(480) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc070msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(1) 
set(gcf,'units','normalized','position',figure_size)%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
ylim([-0.02 0.16])
zlim([-8000 8000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

for j = 1:480
    set(ax, 'XLim', ax.XLim+dt)
    %pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
close(v)

