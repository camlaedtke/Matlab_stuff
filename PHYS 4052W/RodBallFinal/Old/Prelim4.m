% MULTI SPRING WITH FRICTION AND DELAY DIFF EQU

clear all
load('tek0039PrelimParams.mat')
load('SpringSymbolicOdes.mat')

%%
fric = 50;
N = 128;
v0 = N*v_cm;
tspan = [0, 0.005];

v0ball = v0;
y0ball = 0;

Y_0 = zeros(1,2*(N+1));
Y_0(2*N+1) = y0ball;
Y_0(2*N+2) = v0ball;

% wball = 2*(w0*N/pi); % -- 24 masses

wball = (w0*N/pi)/4; % -- 128 masses
% wball = (w0*N/pi)/2; % -- 36 masses
% options = odeset('RelTol',1e-9,'AbsTol',1e-9);
S = ode45(@(t,y)springdde(t,y,w0,wball,N,fric), tspan, Y_0);


t_spr = S.x;

du_left = S.y(1,:);
v_left = S.y(2,:);

du_right = S.y((2*N-1),:)./2;
v_right = S.y(2*N,:);

du_ball = S.y(2*N+1,:);
v_ball = S.y(2*N+2,:);

fprintf('wball = %.4e \n',wball)

%%
t_range = [0.00 0.001];
v_scale = 1*10^6;
offset=0.000005;

close all
font_size = 22;

figure(1)
set(gcf,'units','normalized','position',[0.4 0.5 0.58 0.43],'MenuBar','none');
plot(t_ave(dt), du_ave(dt),':', 'LineWidth', 3, 'Color', 'k','DisplayName','Data')
hold on
plot(t_ave(dt), y_cm,'--','LineWidth', 3 ,'Color', 'black','DisplayName','Vcm')
plot(t_spr+offset, du_right, 'LineWidth', 3,'Color', 'blue','DisplayName','Model')
% plot(t_spr+offset, du_left, 'LineWidth', 3,'Color', 'red')
% plot(t_spr+offset, du_ball, 'LineWidth', 3,'Color', 'g')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim(t_range) 
grid on
set(gca,'FontSize',font_size)
legend

figure(2)
set(gcf,'units','normalized','position',[0.4 0.05 0.58 0.43],'MenuBar','none');
plot(t_ave(dt), vel_ave_filt*v_scale,':', 'LineWidth', 3 ,'Color', 'k',...
    'DisplayName','Data')
hold on
% plot(t_ave(dt), y_v_cm,'--','LineWidth', 3 ,'Color', 'black')
% plot(t_ave(dt), 2*y_v_cm,'--','LineWidth', 3 ,'Color', 'black')
plot(t_spr+offset, v_right, 'LineWidth', 3,'Color', 'blue','DisplayName','Model')
% plot(t_spr+offset, v_left, 'LineWidth', 3,'Color', 'red')
% plot(t_spr+offset, v_ball, 'LineWidth', 3,'Color', 'g')
hold off
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim(t_range) 
grid on
set(gca,'FontSize',font_size)
legend






