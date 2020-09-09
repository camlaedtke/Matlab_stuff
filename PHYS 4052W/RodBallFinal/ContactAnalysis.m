clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')
% load('100MassFirst3ms.mat')
v_no_gravity = diff(du_no_gravity)./diff(t_post_impact);
t_c = 1.2900e-04; % contact time from voltage data\

load('tek0039Sim100M800msecLowTol.mat')

%% Plot All Masses And Contact Force

% Get displacement and velocity vectors for ball, first mass and last mass
t_sim = S.x;
du_m0 = S.y(2*N+1,:);
v_m0 = S.y(2*N+2,:);
du_m1 = S.y(1,:);
v_m1 = S.y(2,:);
du_m100 = S.y((2*N-1),:);
v_m100 = S.y(2*N,:);

% Get information during time where ball & rod are in contact
idx_contact = find((du_m0 - du_m1) > 0);
t_spr_contact = t_sim(idx_contact);
du_ball_contact = du_m0(idx_contact);
du_left_contact = du_m1(idx_contact);
fb = (wball^2)*(du_ball_contact - du_left_contact).^(3/2);

fprintf(' Contact time (model): t_c = %.0f microseconds\n',t_spr_contact(end)*10^6)
fprintf(' Contact time (data): t_c = %.0f microseconds\n\n',t_c*10^6)

% Get vector of N colors ranging from blue to red
r = linspace(0,1,N+2)';
g = zeros(N+2,1);
b = linspace(1,0,N+2)';
rgb = [r,g,b];

% Subtract motion to get in COM frame
du_m1_com = subtractMotion(t_sim, du_m1);
du_m100_com = subtractMotion(t_sim, du_m100);
c0 = v_phase;

%%
% I need to double the value of the x-axis, which keeping the relative
% spacing the same. 

t_sim_norm = t_sim*c0/L;
tc_norm = (t_c/tau)*L;

dt = 1:1400000;
k = (1:N)*2;
t = t_sim_norm(dt);
V = S.y(k,dt);
DU = S.y(k-1,dt);
dv = 1:length(dt)-1;
A = diff(V,1,2)./diff(t_sim(dt));


period_shift = 1.5*1767; % Specify factor to shift x-axis (time axis)
offset = 0.17; 
% 300 msec || 0.17 || 1325
% 200 msec || 0.13 || 882
% 100 msec || 0.07 || 442
t_range = [0 8]; % Specify width of plots in the time dimension


period_shift = period_shift+ offset;
period_shift_msec = (t_range+period_shift+L);
period_shift_msec = period_shift_msec*tau*1000;
t_offset_msec = t_sim;

close all
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.65 0.56 0.3],'MenuBar','none');
plot(t_offset_msec(dt)*2000, V(1,:), 'LineWidth', 2, 'Color', rgb(1,:));
grid on
title('Velocity vs. Time')
xlabel('Time (msec)')
ylabel('Displacement (m)')
xlim(period_shift_msec);
ylim([-0.02 0.12]);
set(gca,'FontSize',18)

t_vel_offset = t_sim_norm - L; 
t_vel_offset = t_vel_offset - period_shift; 

figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.33 0.56 0.3],'MenuBar','none');
plot(t_vel_offset(dt), V(1,:), 'LineWidth', 2, 'Color', 'r');
xlim(t_range)
ylim([-0.02 0.12]);
set(gca,'FontSize',20)
xlabel('Position Along Rod');
ylabel('velocity (m/s)')
title('Velocity vs. Time')
grid on

%%

n_periods = 20;

F(N*n_periods) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('MotionAlongRod100Mfrom400msec', 'MPEG-4');
v.FrameRate = round((3*N)/10);
open(v)

t_range = [0 1]; % Specify width of plots in the time dimension
period_shift = period_shift+ offset;
t_vel_offset = t_sim_norm - L; 
t_accel_offset = t_sim_norm - L;
t_vel_offset = t_vel_offset - period_shift; 
t_accel_offset = t_accel_offset - period_shift;
tau_i = 0;
tau_i = tau_i + period_shift;
dim = [0.46 0.42 0.2 0.1];
% fig_size = [0.42 0.1 0.62 0.62]; (big screen)

close all
figure(1)
set(gcf,'color','w');
fig = gcf;
subplot(2,1,1)
set(gcf,'units','normalized','position',[0.42 0.1 0.62 0.62],'MenuBar','none');
p1 = plot(t_vel_offset(dt), V(1,:), 'LineWidth', 3, 'Color', 'b');
set(p1, 'xdata', t_vel_offset(dt), 'ydata', V(1,:));
ax1 = gca;
ax1.XLim = t_range;
ax1.YLim = [-0.02 0.12];
set(gca,'FontSize',20)
ax1.XLabel.String = 'Position Along Rod';
ax1.YLabel.String = 'velocity (m/s)';
title('Velocity vs. Time')
grid on

subplot(2,1,2)
p2 = plot(t_accel_offset(dv), A(1,:), 'LineWidth', 3, 'Color', 'b');
set(p2, 'xdata', t_accel_offset(dv), 'ydata', A(1,:));
ax2 = gca;
ax2.XLim = t_range;
ax2.YLim = [-2000 2000];
% 300 msec || 2000
% 200 msec || 2100 
% 100 msec || 2200 
set(gca,'FontSize',20)
ax2.XLabel.String = 'Position Along Rod';
ax2.YLabel.String = 'Acceleration (m/s^2)';
str = ['Time: ',num2str(round(tau*tau_i*1000,1)),' msec'];
text = annotation(fig,'textbox',dim,'String',str,'FitBoxToText','on');
text.FontSize = 28;
title('Acceleration vs. Time')
grid on

for i = 1:n_periods*N
    if mod(i,N) == 0 && i ~= 0
        tau_i = tau_i + 1;
        text.String = ['Time:  ', num2str(round(tau*tau_i*1000,1)),'  msec'];
        t_vel_offset = t_vel_offset - 1; % try -2
        t_accel_offset = t_accel_offset - 1; % try -2
    end
    n = mod(i,N)+1;
    set(p1, 'xdata', t_vel_offset(dt), 'ydata', V(n,:));
    set(p2, 'xdata', t_accel_offset(dv), 'ydata', A(n,:));  
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
    pause(0.03)
end
close(v)

%% 2-D Plots

t_range = [0.00 2*tau]; % Specify width of plots in the time dimension
xlim_shift = 0.0; % Specify factor to shift x-axis (time axis)
font_size = 18;

% create time offset
dt = tau/2;
% dt = 2.0526e-04;
dt_n = zeros(1,100);
p = linspace(0,1,100);
for j = 1:100
    dt_n(j) = p(j)*dt; %*v_cm*tau/2;
end

% Plot of displacement vs time for all 100 masses simultaniously
close all
figure(1)
set(gcf,'units','normalized','position',[0.42 0.65 0.56 0.3],'MenuBar','none');
% plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
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
xlim(t_range+xlim_shift) 

% % Plot of velocity vs time for all 100 masses simultaniously
figure(2)
set(gcf,'units','normalized','position',[0.42 0.34 0.56 0.3],'MenuBar','none');
% plot(t_post_impact(dt), v_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
hold on
plot(t_sim(idx_contact), 2*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot(t_sim,S.y(i*2,:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim(t_range+xlim_shift) 
set(gca,'FontSize',font_size)

% Plot of acceleration vs time for all 100 masses simultaniously
figure(3)
set(gcf,'units','normalized','position',[0.42 0.03 0.56 0.3],'MenuBar','none');
% plot(t_post_impact(dt), v_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
hold on
plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot(t_sim(1:end-1),diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Acceleration (m/s^2)')
xlim(t_range+xlim_shift) 
set(gca,'FontSize',font_size)



%% 3-D Plots

t_range = [0.00 tau];
xlim_shift = 0;

close all
font_size = 18;
figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44])%,'MenuBar','none');
% plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
hold on
% plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot3(t_sim,S.y((i*2-1),:),S.y((i*2),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
zlabel('Velocity (m/s)')
grid on
view(290,15)
set(gca,'FontSize',font_size)
xlim(t_range+xlim_shift) 


% Define video object 
% F(480) = struct('cdata',[],'colormap',[]); 
% v = VideoWriter('BallRodContactvid', 'MPEG-4');
% v.FrameRate = 24;
% open(v)


figure(2)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.8])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot3(t_sim(1:end-1),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
% ylim([-0.01 0.1])
% zlim([-2500 2500])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-50,15)

% dt = tau/240;
% for j = 1:480
%     set(ax, 'XLim', ax.XLim+dt)
%     % psi = psi-dpsi;
%     % view(psi,15)
%     pause(2);
%     F(i) = getframe(gcf);
%     writeVideo(v, F(i))
% end
% close(v)


%% Plot Motion of All Masses During Contact Time

t_range = [0.00 0.4*10^-5];
xlim_shift = 0;
close all
font_size = 18;
figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44])%,'MenuBar','none');
plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
hold on
plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot(t_sim(idx_contact),S.y((i*2-1),idx_contact),'LineWidth',1,'Color',rgb(i,:))
end
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim(t_range+xlim_shift) 
% ylim([-1*10^-5 6*10^-5])
grid on
set(gca,'FontSize',font_size)
% legend

% figure(2)
% set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44])%,'MenuBar','none');
% plot(t_post_impact(dt), v_no_gravity(dt),':', 'LineWidth', 2 ,'Color', 'k')
% hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
% for i=1:100
%     plot(t_sim(idx_contact),S.y(i*2,idx_contact),'LineWidth',1,'Color',rgb(i,:))
% end
% hold off
% title('Velocity vs. Time')
% xlabel('Time (seconds)')
% ylabel('Velocity (m/s)')
% xlim(t_range+xlim_shift) 
% % ylim([-0.25 0.25])
% grid on
% set(gca,'FontSize',font_size)

figure(2)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44])%,'MenuBar','none');
plot(t_post_impact(dt), v_no_gravity(dt),':', 'LineWidth', 2 ,'Color', 'k')
hold on
plot(t_sim(idx_contact), 0.05*fb,'LineWidth',4,'Color', 'g')
for i=1:100
    dv = diff(S.y(i*2,:))./diff(t_sim);
    plot(t_sim(idx_contact),dv(idx_contact),'LineWidth',1,'Color',rgb(i,:))
end
hold off
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Acceleration (m/s)')
xlim(t_range+xlim_shift) 
% ylim([-0.25 0.25])
grid on
set(gca,'FontSize',font_size)

%% ------ Plots with 'fake' displacement added -------

dx = v_cm*tau/2;
% dx = 10e-06;
dx_n = zeros(1,100);

p = linspace(0,1,100);

X = zeros(200,length(du_m1));
for i = 1:100
    dx_n(i) = p(i)*dx;
    X((i*2-1),:) = S.y((i*2-1),:) + dx_n(i);
    X(i*2,:) = S.y(i*2,:);
end

% Define time shift per frame
% Define video object 
% F(200) = struct('cdata',[],'colormap',[]); 
% v = VideoWriter('BallRodVid1', 'MPEG-4');
% v.FrameRate = 15;
% open(v)

t_range = [0.00 tau];
xlim_shift = 0.79; 
font_size = 15;

close all
figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44])%,'MenuBar','none');
hold on
plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot(t_sim,X((i*2-1),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize',font_size) 
ax1 = gca;
ax1.XLim = t_range+xlim_shift;

figure(2)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44])%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim,X((i*2-1),:),X((i*2),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
zlabel('Velocity (m/s)')
set(gca,'FontSize',font_size)
ax2 = gca;
ax2.XLim = t_range+xlim_shift;
view(-10,15)

% df = 0.1*tau;
% for j = 1:200
%     set(ax1, 'XLim', ax1.XLim+df)
%     set(ax2, 'XLim', ax2.XLim+df)
%     %view(290,15)
%     pause(3);
%     F(i) = getframe(gcf);
%     writeVideo(v, F(i))
% end
% close(v)

%% ------ Plots with time offset added --------

% create time offset
dt = tau/2;
% dt = 2.0526e-04;
dt_n = zeros(1,100);
p = linspace(0,1,100);
for j = 1:100
    dt_n(j) = p(j)*dt; %*v_cm*tau/2;
end

t_range = [0.00 tau];
xlim_shift = -tau; 
font_size = 15;

close all
% figure(1) 
% set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44])%,'MenuBar','none');
% % plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 2,'Color', 'k')
% hold on
% % plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
% for i=1:100
%     plot(t_sim+dt_n(i),S.y((i*2-1),:),'LineWidth',1,'Color',rgb(i,:))
% end
% hold off
% grid on
% title('Displacement vs. Time')
% xlabel('Time (seconds)')
% ylabel('Displacement (m)')
% set(gca,'FontSize',font_size)
% xlim(t_range+xlim_shift) 

% figure(2)
% set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44])%,'MenuBar','none');
% hold on
% % plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
% for i=1:100
%     plot(t_sim+dt_n(i),S.y((i*2),:),'LineWidth',1,'Color',rgb(i,:))
% end
% hold off
% grid on
% title('Velocity vs. Time')
% xlabel('Time (seconds)')
% ylabel('Displacement (m)')
% set(gca,'FontSize',font_size)
% xlim(t_range+xlim_shift) 

% figure(3) 
% set(gcf,'units','normalized','position',[0.42 0.05 0.28 0.44])%,'MenuBar','none');
% hold on
% for i=1:100
%     plot3(t_sim+dt_n(i),S.y((i*2-1),:),S.y((i*2),:),'LineWidth',1,'Color',rgb(i,:))
% end
% hold off
% grid on
% title('Displacement vs. Velocity vs. Time')
% xlabel('Time (seconds)')
% ylabel('Displacement (m)')
% zlabel('Velocity (m/s)')
% set(gca,'FontSize',font_size)
% xlim(t_range+xlim_shift) 
% view(20,15)

% Define video object 

F(720) = struct('cdata',[],'colormap',[]); 
v = VideoWriter('BallRodFollowTc080msec', 'MPEG-4');
v.FrameRate = 24;
open(v)

figure(4) 
set(gcf,'units','normalized','position',[0.7 0.05 0.56 0.8])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot3(t_sim(1:end-1)+dt_n(i),S.y(i*2,1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim+dt_n(i)),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
% xlim(t_range+xlim_shift) 
ylim([-0.02 0.16])
zlim([-6000 6000])
ax = gca;
ax.XLim = t_range+xlim_shift;
view(-85,15)

% Animation parameters
% psi = -15;
dt = tau/240;

% Animate first period or so
for j = 1:720
    set(ax, 'XLim', ax.XLim+dt)
    % psi = psi+dpsi;
    % view(psi,15)
    pause(2);
    F(i) = getframe(gcf);
    writeVideo(v, F(i))
end
% 
% % Jump ahead to 0.1 sec
% set(ax, 'XLim', ax.XLim+0.1)
% for j = 1:480
%     set(ax, 'XLim', ax.XLim+dt)
%     % psi = psi+dpsi;
%     % view(psi,15)
%     pause(2);
%     F(i) = getframe(gcf);
%     writeVideo(v, F(i))
% end
% 
% % Jump ahead to 0.2 sec
% set(ax, 'XLim', ax.XLim+0.1)
% for j = 1:480
%     set(ax, 'XLim', ax.XLim+dt)
%     % psi = psi+dpsi;
%     % view(psi,15)
%     pause(2);
%     F(i) = getframe(gcf);
%     writeVideo(v, F(i))
% end
% 
close(v)



%% 2-D strain wave propagation

close all
figure(1)
set(gcf,'units','normalized','position',[0.42 0.02 0.56 0.88])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot(t_sim(1:end-1)+dt_n(i),...
        diff(S.y(i*2,:))./diff(t_sim+dt_n(i)),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim(t_range+xlim_shift) 
set(gca,'FontSize',font_size)


%% 3-D strain wave propagation

close all
figure(1)
set(gcf,'units','normalized','position',[0.42 0.02 0.56 0.88])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot3(t_sim(1:end-1)+dt_n(i),S.y((i*2),1:end-1),...
        diff(S.y(i*2,:))./diff(t_sim+dt_n(i)),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
xlim(t_range+xlim_shift) 
set(gca,'FontSize',font_size)
view(20,15)


%% ------ Plots with displacement and time offset -------

dx = v_cm*tau/2;
dt = tau/2;
% dx = 10e-06;
% dt = 0.1e-04;

dx_n = zeros(1,100);
dt_n = zeros(1,100);
p = linspace(0,1,100);

X = zeros(200,length(du_m1));
for i = 1:100
    dx_n(i) = p(i)*dx; % create position offset
    dt_n(i) = p(i)*dt; % create time offset
    X((i*2-1),:) = S.y((i*2-1),:) + dx_n(i); 
    X(i*2,:) = S.y(i*2,:);
end

t_range = [0.00 tau];
xlim_shift = 0.79; 
font_size = 15;

close all

figure(1) 
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44])%,'MenuBar','none');
hold on
plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot(t_sim+dt_n(i),X((i*2-1),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize',font_size)
xlim(t_range+xlim_shift) 

% Plot of velocity vs time
figure(2) 
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44])%,'MenuBar','none');
hold on
plot(t_sim(idx_contact), 0.4*fb*10^-10,'LineWidth', 4,'Color', 'g')
for i=1:100
    plot(t_sim+dt_n(i),X((i*2),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize',font_size)
xlim(t_range+xlim_shift) 

figure(3) 
set(gcf,'units','normalized','position',[0.42 0.03 0.28 0.44])%,'MenuBar','none');
hold on
for i=1:100
    plot3(t_sim+dt_n(i),X((i*2-1),:),X((i*2),:),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Displacement vs. Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
zlabel('Velocity (m/s)')
set(gca,'FontSize',font_size)
xlim(t_range+xlim_shift) 
view(20,15)

figure(4) 
set(gcf,'units','normalized','position',[0.7 0.03 0.28 0.44])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot3(t_sim(1:end-1)+dt_n(i),X(i*2,1:end-1),...
        diff(X(i*2,:))./diff(t_sim+dt_n(i)),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Velocity vs. Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
zlabel('Acceleration (m/s^2)')
set(gca,'FontSize',font_size)
xlim(t_range+xlim_shift) 
view(20,15)

%% ------ Plots with velocity offset -------

dx = 0;
dt = 0;
dv = 0;
% dx = 10e-06;
% dt = 0.1e-04;

dx_n = zeros(1,100);
dt_n = zeros(1,100);
dv_n = zeros(1,100);
p = linspace(0,1,100);

X = zeros(200,1352661);
for i = 1:100
    dx_n(i) = p(i)*dx; % create position offset
    dt_n(i) = p(i)*dt; % create time offset
    dv_n(i) = p(i)*dv;
    X((i*2-1),:) = S.y((i*2-1),:) + dx_n(i); 
    X(i*2,:) = S.y(i*2,:) + dv_n(i);
end

t_range = [0.00 tau];
xlim_shift = tau + 3*tau/4; 
font_size = 15;

close all
figure(2)
set(gcf,'units','normalized','position',[0.42 0.02 0.56 0.88])%,'MenuBar','none');
hold on
% plot(t_sim(idx_contact), 1.5*fb*10^-6,'LineWidth',4,'Color', 'g')
for i=1:100
    plot(t_sim(1:end-1)+dt_n(i),...
        diff(X(i*2,:))./diff(t_sim+dt_n(i)),'LineWidth',1,'Color',rgb(i,:))
end
hold off
grid on
title('Acceleration vs. Time')
xlabel('Time (seconds)')
ylabel('Acceleration (m/s^2)')
xlim(t_range+xlim_shift) 
set(gca,'FontSize',font_size)








