% Prelim 5
clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')
dt = 1:490000;
v_no_gravity = diff(du_no_gravity)./diff(t_post_impact);
t_c = 1.2900e-04; % contact time from voltage data

%%
fric = 8;
N = 100; 

v0 = N*v_cm/2; % Changes particle velocity and CM velocity
lambda = 800000; % Changes characteristic frequency
w0_n = (w0*N)/pi;
kball = 4.84*10^9;
wball = sqrt(kball/m_ball);

tspan = [0, 0.0005];

v0ball = 1.023*v0; % tweak initial velocity 
y0ball = 0;

Y_0 = zeros(1,2*(N+1));
Y_0(2*N+1) = y0ball;
Y_0(2*N+2) = v0ball;

tic
options = odeset('stats','on','RelTol',1e-13,'AbsTol',1e-13); %'OutputFcn',@odeplot,
S = ode45(@(t,y)springdde(t,y,w0_n,wball,lambda,N,fric), tspan, Y_0,options);
toc

S = odextend(S,@(t,y)springdde(t,y,w0_n,wball,lambda,N,fric),0.001);



%%

t_spr = S.x;
du_left = S.y(1,:);
v_left = S.y(2,:);
du_middle = S.y(N-1,:);
v_middle = S.y(N,:);
du_right = S.y((2*N-1),:);
v_right = S.y(2*N,:);
du_ball = S.y(2*N+1,:);
v_ball = S.y(2*N+2,:);

% load('tek0039Simulation72Springs.mat')
idx_contact = find((du_ball - du_left) > 0);
t_spr_contact = t_spr(idx_contact);
du_ball_contact = du_ball(idx_contact);
du_left_contact = du_left(idx_contact);
fb = (wball^2)*(du_ball_contact - du_left_contact).^(3/2);
%
% INTEGRATE FB, use m_ball,  m_rod, v_cm to get v_0?
fprintf(' Contact time (model): t_c = %.0f microseconds\n',t_spr_contact(end)*10^6)
fprintf(' Contact time (data): t_c = %.0f microseconds\n\n',t_c*10^6)

t_range = [0.00 0.001];
xlim_shift = 0;

close all
font_size = 18;
figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44],'MenuBar','none');
plot(t_spr, du_right, 'LineWidth', 3,'Color', 'r','DisplayName','right')
hold on
plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 3,'Color', 'k','DisplayName','Data')
% plot(t_spr, du_left, 'LineWidth', 3,'Color', 'b','DisplayName','left')
% plot(t_spr, du_middle, 'LineWidth', 3,'Color', 'c','DisplayName','middle')
% plot(t_spr(idx_contact), du_ball(idx_contact)*0.2, 'LineWidth', 3,'Color', 'g','DisplayName','ball')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim(t_range+xlim_shift) 
grid on
set(gca,'FontSize',font_size)
legend

figure(2)
set(gcf,'units','normalized','position',[0.42 0.05 0.28 0.44],'MenuBar','none');
plot(t_spr, v_right,'LineWidth',3 ,'Color', 'r','DisplayName','right')
hold on
plot(t_post_impact(dt), v_no_gravity(dt),':', 'LineWidth', 3 ,'Color', 'k','DisplayName','Data')
% plot(t_spr, v_left, 'LineWidth', 3,'Color', 'b','DisplayName','left')
% plot(t_spr, v_middle, 'LineWidth', 1,'Color', 'c','DisplayName','middle')
% plot(t_spr(idx_contact), v_ball(idx_contact)*0.1, 'LineWidth', 3,'Color', 'g','DisplayName','ball')
hold off
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim(t_range+xlim_shift) 
grid on
set(gca,'FontSize',font_size)
legend

figure(3)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44],'MenuBar','none');
plot(t_spr, du_right, 'LineWidth', 3,'Color', 'r','DisplayName','right')
hold on
plot(t_post_impact(dt), du_no_gravity(dt),':', 'LineWidth', 3, 'Color', 'k','DisplayName','Data')
% plot(t_spr+offset, du_left, 'LineWidth', 3,'Color', 'b','DisplayName','left')
% plot(t_spr+offset, du_middle, 'LineWidth', 1,'Color', 'c','DisplayName','middle')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
% xlim([0 0.01])
grid on
set(gca,'FontSize',font_size)
legend

figure(4)
set(gcf,'units','normalized','position',[0.7 0.05 0.28 0.44],'MenuBar','none');
patchline(t_post_impact(dt), v_no_gravity(dt),'edgecolor','k','linewidth',0.05,'edgealpha',0.05)
hold on
patchline(t_spr, v_right,'edgecolor','r','linewidth',0.05,'edgealpha',0.05)
% patchline(t_spr, v_left,'edgecolor','b','linewidth',0.01,'edgealpha',0.03)
% patchline(t_spr, v_middle,'edgecolor','c','linewidth',0.01,'edgealpha',0.01)
hold off
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim([0 1]) 
grid on
set(gca,'FontSize',font_size)

