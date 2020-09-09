
clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')

dt = 1:490000;
v_no_gravity = diff(du_no_gravity)./diff(t_post_impact);
t_c = 1.2900e-04; % contact time from voltage data

% 
% Get damping coefficient from a logrithmic fit to the data of 
% displacement vs. time with pendulum subtracted
% 1) Import tek0039PendulumRemoved.mat into smoothing script
% 2) Grab displacment data of free-end
% 3) Put it through a lowpass filter 
%      - If it is low enough, I might get the curve for a damped 
%        harmonic oscillator.
% 4) Use findpeaks to get the amplitude of the first peak, and second peak
% 5) According to bookmarked page, it looks something like that
%      - See if 1/n ln|x1/x_n+1| changes when you compare the first peak
%        to the second, to the third, etc. 
% If that doesn't
% - No lowpass filter. Just use findpeaks on the displacement data
% - Save the findpeaks dataset to a matfile. 
% Make the following algorithm
% For Try first 300 msec of data, get 
%   

beta = -3.41;
w1 = sqrt(w0^2-beta^2);
fprintf('w0: %.6f kHz |  w1: %.6f kHz \n', w0*10^-3, w1*10^-3)

%%
fric = 4;

N = 64; 

v0 = N*v_cm/2; % Changes particle velocity and CM velocity

lambda = 0;
w0_n = (w0*N)/pi; 
% w_damp = sqrt(w0_n^2 - beta^2);
% kball = 4.84*10^9;
kball = 5*10^9;
wball = sqrt(kball/m_ball);

tspan = [0,0.5];

v0ball = 1.01*v0; % tweak initial velocity 
y0ball = 0;
Y_0 = zeros(1,2*(N+1));
Y_0(2*N+1) = y0ball;
Y_0(2*N+2) = v0ball;

fprintf(['\n\n ---------- Input Parameters ---------- \n'...
        ' w0_n   = %.0f rad/sec ||| f0_n = %.0f Hz \n'...
        ' w_ball = %.0f rad/sec  ||| v0   = %.0f v_cm \n'...
     '       N = %.0f masses     ||| N-1  = %.0f springs \n'],...
            w0_n,w0_n/(2*pi),wball,v0ball/v_cm,N,(N-1))

tic
% options = odeset('stats','on','RelTol',1e-12,'AbsTol',1e-12); 
options = odeset('stats','on','RelTol',1e-10,'AbsTol',1e-10); 
S = ode45(@(t,y)springdde(t,y,w0_n,wball,lambda,N,fric), tspan, Y_0,options);
toc

%
% save('tek0039Sim256M5msecLowTol.mat','S','N','w0_n','v0','v0ball',...
%              'wball','lambda','fric','-v7.3')
%
%
%
% load('tek0039Simulation100SpringsLowTol.mat')
% tic 
% S = odextend(S,@(t,y)springdde(t,y,w0_n,wball,lambda,N,fric),0.9);
% toc
% save('tek0039Simulation100SpringsLowTol.mat','S','N','w0_n','v0','v0ball',...
%                                  'wball','lambda','fric','-v7.3')


t_spr = S.x;
du_left = S.y(1,:);
v_left = S.y(2,:);
du_middle = S.y(N-1,:);
v_middle = S.y(N,:);
du_right = S.y((2*N-1),:);
v_right = S.y(2*N,:);
du_ball = S.y(2*N+1,:);
v_ball = S.y(2*N+2,:);
%
% load('tek0039Simulation100Springs.mat')

idx_contact = find((du_ball - du_left) > 0);
t_spr_contact = t_spr(idx_contact);
du_ball_contact = du_ball(idx_contact);
du_left_contact = du_left(idx_contact);
fb = (wball^2)*(du_ball_contact - du_left_contact).^(3/2);
%
% INTEGRATE FB, use m_ball,  m_rod, v_cm to get v_0?
fprintf('Contact time (model): t_c = %.0f microseconds\n',t_spr_contact(end)*10^6)
fprintf('Contact time (data): t_c = %.0f microseconds\n\n',t_c*10^6)

t_range = [0.00 4*tau];
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
% plot(t_spr(idx_contact), v_ball(idx_contact)*0.1, 
% 'LineWidth', 3,'Color', 'g','DisplayName','ball')
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

y_vcm = v_cm*ones(length(t_post_impact(dt)),1);
figure(4)
set(gcf,'units','normalized','position',[0.7 0.05 0.28 0.44],'MenuBar','none');
patchline(t_post_impact(dt), v_no_gravity(dt),'edgecolor','k','linewidth',0.05,'edgealpha',0.05)
patchline(t_post_impact(dt), y_vcm,'edgecolor','k','linewidth',1,'edgealpha',1)
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


% Propagated Stress / Force Pulse
t_scale = 10^6;
% figure(5)
% set(gcf,'units','normalized','position',[0.3 0.3 0.28 0.44]);
% hold on
% plot(t_scale*t_spr_contact, fb,...
%     'LineWidth',3 ,'Color', 'c','DisplayName','Contact Force (Model Data)')
% plot( t_scale*t_spr_contact, du_left_contact*10^9, ...
%      ':','LineWidth', 3,'Color', 'b','DisplayName','Rod')
% plot( t_scale*t_spr_contact, du_ball_contact*10^9, ...
%     ':','LineWidth', 3,'Color', 'k','DisplayName','Ball')
% hold off
% title('Contact Force vs. Time')
% xlabel('Time ($\mu s$)','Interpreter','latex')
% ylabel('Contact Force (N)','Interpreter','latex')
% grid on
% set(gca,'FontSize',font_size)
% legend('Location','northeast','Interpreter','latex')
 
%% Interpolate model data 

dt = 1:40000;

t_spr_interp = linspace(0,0.1,1400000);
du_right_interp = interp1(t_spr,du_right,t_spr_interp);
v_right_interp = interp1(t_spr,v_right,t_spr_interp);
v_middle_interp = interp1(t_spr,v_middle,t_spr_interp);
font_size = 18;

close all

[spr_velFT_right,  spr_f_right, Fs] =  MyFFT(v_right_interp,t_spr_interp); 
[velFT,f] = MyFFT(v_no_gravity(dt), t_post_impact(dt)); 

figure(5)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44])%,'MenuBar','none');
plot(f, velFT,'Color','k','linewidth',4)
hold on
plot(spr_f_right, spr_velFT_right,'Color','r','linewidth',2)
hold off
title('Fourier Transform of the Data and Model (Velocity)')
xlabel('f (Hz)')
ylabel('Amplitude')
xlim([0, 350000])
grid on
legend('Model', 'Data')
set(gca,'FontSize',font_size,'YScale', 'log')
fprintf('Sampling Freq: Fs = %.10e \n', Fs)

%%
close all
[~,~, spr_Fs] = MyFFT(v_right_interp,t_spr_interp); 
[~,~, Fs] = MyFFT(v_no_gravity(dt), t_post_impact(dt)); 


L = size(v_no_gravity(dt),1);
res = Fs/L;
% Spectrogram - velocity data
figure(8)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44]);
[sp,fp,tp] = pspectrum(v_no_gravity(dt),Fs,'spectrogram','FrequencyResolution',200*res);
mesh(tp,fp,pow2db(sp))
% mesh(tp,fp,sp)
% waterfall(tp,fp,sp)
view(95,15)
title('Spectrogram (Velocity Data Pendulum Subtracted)')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([300 20000])
% zlim([0 1.2*10^-3])
zlim([-140 -20])
zlabel('Power')
set(gca,'FontSize',15)

% close all
L = size(v_right_interp,2);
res = spr_Fs/L;
figure(9)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44]);
[sp,fp,tp] = pspectrum(v_right_interp,spr_Fs,'spectrogram','FrequencyResolution',200*res);
mesh(tp,fp,pow2db(sp))
% mesh(tp,fp,sp)
% waterfall(tp,fp,sp)
view(95,15)
title('Spectrogram (100 Mass with Friction and Hertz)')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([300 20000])
% zlim([0 10^-3])
zlim([-200 -20])
zlabel('Power')
set(gca,'FontSize',15)

%%
close all

figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44]);
[imf,~,~] = emd(v_no_gravity(dt),'Interpolation','pchip');
[~,~,t,~,~] = hht(imf,Fs);
[p,q] = ndgrid(t,1:size(imf,2));
plot3(p,q,imf,'linewidth',3)
grid on
xlabel('Time Values')
ylabel('Mode Number')
zlabel('Mode Amplitude')
set(gca,'FontSize',18)

figure(2)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
pspectrum(v_no_gravity(dt),Fs,'persistence','FrequencyLimits',[0 20000],'TimeResolution',0.01)
set(gca,'FontSize',18)

% You only see cool stuff when signal is super filtered (velocity)
figure(3)
set(gcf,'units','normalized','position',[0.42 0.05 0.28 0.44]);
[imf,~,~] = emd(v_right_interp,'Interpolation','pchip');
[~,~,t,~,~] = hht(imf,spr_Fs);
[p,q] = ndgrid(t,1:size(imf,2));
plot3(p,q,imf,'linewidth',3)
grid on
xlabel('Time Values')
ylabel('Mode Number')
zlabel('Mode Amplitude')
set(gca,'FontSize',18)

figure(4)
set(gcf,'units','normalized','position',[0.7 0.05 0.28 0.44]);
pspectrum(v_right_interp,spr_Fs,'persistence','FrequencyLimits',[0 20000],'TimeResolution',0.01)
set(gca,'FontSize',18)

%%
save('tek0039Sim200M50msecLowTol.mat','fric','v0ball','kball','wball',...
    'dt','t_spr',...
    'idx_contact','t_spr_contact','fb',...
    'du_right','du_left','du_middle','du_ball',...
    'v_right','v_left','v_middle','v_ball')









