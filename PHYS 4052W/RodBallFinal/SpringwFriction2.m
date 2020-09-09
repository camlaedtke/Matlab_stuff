clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')

t_data = t_post_impact; 
u_data = du_no_gravity; 
v_data = v_no_gravity;
c = v_phase;

%% Main model parameters
N = 100; 
beta = 3.334; % from fit to envelope peaks
tspan = [0, 0.3];
kball = 1.5*10^10;
wball = sqrt(kball/m_ball);
% v0ball = N*v_cm/2; % wrong v0, but gives best results
v0ball = (v_cm*(m_ball+m_rod))./(2*m_ball); % correct v0
y0ball = 0;
w1 = sqrt(w0^2 - beta^2);
w0_n = (w1*N)/pi;

Y_0 = zeros(1,2*(N+1));
Y_0(2*N+1) = y0ball;
Y_0(2*N+2) = v0ball;

options = odeset('stats','on','RelTol',1e-11,'AbsTol',1e-11); 
S = ode45(@(t,y)springdde(t,y,w0_n,wball,N,beta), tspan, Y_0, options);
% Extend ode45 solution. Load/Save ode45 solution
% load('tek0039Sim100MFullSecFixedV0.mat')
% options = odeset('stats','on','RelTol',1e-9,'AbsTol',1e-10); 
% S = odextend(S,@(t,y)springdde(t,y,w0_n,wball,N,beta),1.5,[],options);
% save('tek0039Sim128M800msecFixedV0.mat','S','N','w0_n','v0ball','wball','kball','beta','-v7.3')
%
t_sim = S.x';
u_left = S.y(1,:)';
v_left = S.y(2,:)';
u_right = S.y((2*N-1),:)';
v_right = S.y(2*N,:)';
u_ball = S.y(2*N+1,:)';
v_ball = S.y(2*N+2,:)';
[idx_contact, t_contact, fb] = contactdata(u_ball, u_left, t_sim, wball);

SCALE = 9.60089; 
dx = N/L_rod; k = (1:N)*2;

U = S.y(k-1,:);
DU = diff(U); % differentiated in row direction (length)
DUDX = DU./dx; DUDX(N,:) = zeros(1,length(t_sim)); % zero pad last row

% Display results
plot_params = struct('linewidth',0.05,'edgealpha',0.05);
plot_params.tspan = [0, 1];
plot_params.font_size = 20;
close all
displacementresponse(t_sim, u_right*SCALE, v_right*SCALE,...
    t_data, u_data, v_data, plot_params)

%% Get Strain
M = m_rod/N;
F_c = -M*diff(sgolayfilt(v_data,3,21));

% Stress
stress_sim_left = Y*DUDX(1,:);
stress_sim_right = Y*DUDX(N-1,:);
stress_data = F_c/A_rod; stress_data(end+1)=0;

% Strain
strain_sim_left = DUDX(1,:);
strain_sim_right = DUDX(N-1,:);
strain_data = stress_data*(1/Y);

close all
tiledlayout(3,1);
fig_size = [0.42 0.05 0.58 0.9];
t_range = [0 1];
t_shift = 0;
font_size = 22;

nexttile
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot(t_sim*10^3, u_right*SCALE,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, u_data,':','Color','k','LineWidth',3)
hold off
grid on
title('Displacement (Scaled)')
xlabel('Time [msec]','interpreter','latex')
ylabel('Displacement [meters]','interpreter','latex')
xlim(t_range+t_shift)
% ylim([-0.015 0.1])
set(gca,'FontSize',font_size)

nexttile
plot(t_sim*10^3, v_right*SCALE,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
title('Particle Velocity at Free End (Scaled)')
xlabel('Time [msec]','interpreter','latex')
ylabel('Velocity [m/s]','interpreter','latex')
xlim(t_range+t_shift)
% ylim([-0.015 0.1])
set(gca,'FontSize',font_size)

nexttile
plot(t_sim*10^3, strain_sim_right,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, strain_data,':','Color','k','LineWidth',3)
hold off
grid on
title('Strain (Not Scaled)')
xlabel('Time [msec]','interpreter','latex')
ylabel('Strain','interpreter','latex')
xlim(t_range+t_shift)
% ylim([-0.015 0.1])
set(gca,'FontSize',font_size)


%% Plot Comparason of strain/velocity waveforms

close all
plot_params.fig_size = [0.42 0.5 0.56 0.44];
plot_params.data_linewidth = 0.3;
plot_params.model_linewidth = 0.3;
plot_params.t_range = [0 1];
plot_params.xlim_shift = 0;
plot_params.font_size = 32;
% velocityplots(t_sim, v_right, t_data, v_data, plot_params)
plot_params.ylim = [-0.2e-11 0.2e-11];
strainplots(t_sim, strain_sim_right,t_data, strain_data,plot_params)

%% Strain/Stress/Velocity Plots
plot_params.t_range = [0 1];
plot_params.xlim_shift = 0;
plot_params.data_linewidth = 0.3;
plot_params.model_linewidth = 0.3;
plot_params.font_size = 28;
plot_params.line_type = ':';

close all
tiledlayout(3,1);
plot_params.fig_size = [0.42 0.05 0.56 0.9];
plot_params.ylim = [-0.3 0.3];
nexttile
stressplots(t_sim, stress_sim_right,t_data, stress_data,plot_params)
nexttile
velocityplots(t_sim, v_right*SCALE, t_data, v_data, plot_params)
plot_params.data_linewidth = 3;
plot_params.model_linewidth = 3;
nexttile
displacementplots(t_sim, u_right*SCALE, t_data, u_data, plot_params)

%% Load Unaveraged Raw Data for FFT's
load('tek0039RawDisplacementAndVelocity.mat')
t_data = t;
u_data = du;
v_data = dudt; v_data(end+1) = 0;

%% Check filtering of raw data
v_data_filt = sgolayfilt(v_data,5,331);
F_c = -m_rod*diff(sgolayfilt(v_data_filt,2,31));
F_c_raw = -m_rod*diff(v_data_filt);
strain_data = F_c./A_rod; strain_data(end+1) = 0;
strain_data_raw = F_c_raw./A_rod; strain_data_raw(end+1) = 0;

%%
fig_size = [0.42 0.1 0.58 0.8];
t_range = [0.0049 0.0054];
t_shift = 0.0;
font_size = 22;

close all
tiledlayout(2,1);
nexttile
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot(t_data, v_data,':','Color','k','LineWidth',1)
hold on
plot(t_data, v_data_filt,'Color','g','LineWidth',2)
hold off
grid on
title('Velocity')
xlabel('Time [sec]','interpreter','latex')
ylabel('Velocity [m/s]','interpreter','latex')
xlim(t_range+t_shift)
ylim([-0.02 0.13])
set(gca,'FontSize',font_size)

nexttile
% plot(t_data, strain_data_raw,':','Color','k','LineWidth',0.3)
% hold on
plot(t_data, strain_data,'Color','g','LineWidth',3)
hold on
plot(t_sim+0.00487, stress_sim_right*7,'Color','r','LineWidth',3)
hold off
grid on
title('Stress (Scaled)')
xlabel('Time [sec]','interpreter','latex')
ylabel('Stress [$N/m^{2}$]','interpreter','latex')
xlim(t_range+t_shift)
% ylim([-1.4 1.4])
set(gca,'FontSize',font_size)



%% Get frequency
% Interpolate model data 
t_sim_interp = linspace(0,0.3,length(t_sim)); % get even spacing btwn points
v_right_interp = interp1(t_sim, v_right*SCALE, t_sim_interp); 

[sim_FFT,  sim_freqs, sim_fs] =  MyFFT(v_right_interp,t_sim_interp); 
[data_FFT,data_freqs, data_fs] = MyFFT(v_data_filt, t_data); % use filtered 
sim_freqs_khz = sim_freqs./1000; % Get in kHz for better plotting
data_freqs_khz = data_freqs./1000;
sim_freqs_norm = 2*pi*sim_freqs./w0;
data_freqs_norm = 2*pi*data_freqs./w0;

%% Compare frequency spectra
plot_params.freq_range = [0, 300];
plot_params.freq_shift = 0;
plot_params.ylim = [1*10^-10 2*10^-2];
plot_params.data_linewidth = 2;
plot_params.model_linewidth = 3;
plot_params.yscale = 'log'; % 'linear'
plot_params.xlabel = 'Frequency [kHz]';
plot_params.title = 'Frequency Spectra';
plot_params.font_size = 22;

close all
plot_params.fig_size = [0.2 0.5 0.8 0.5];
frequencyplots(data_freqs_khz, data_FFT,sim_freqs_khz,sim_FFT,plot_params)

%% Compare "normalized" frequency spectra
close all
plot_params.freq_range = [0, 8];
plot_params.freq_shift = 0;
plot_params.ylim = [2*10^-8 2*10^-2];
plot_params.xlabel = '$\omega_0$';
plot_params.title = 'FFT of Velocity at Free-End';
frequencyplots(data_freqs_norm, data_FFT,sim_freqs_norm,sim_FFT,plot_params)

%% Persistence Spectra: Plots all calculated FFT's over given time period

v_data_filt = sgolayfilt(v_data,5,331);
font_size = 22;
f_min = 0;
f_lim = 500000;
t_res = 0.02;
close all

tiledlayout(3,1);
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.58 0.9],'MenuBar','none');
nexttile
pspectrum(v_data, data_fs,'persistence',...
    'FrequencyLimits',[f_min f_lim],'TimeResolution',t_res)
title('Unfiltered Raw Data')
set(gca,'FontSize',font_size)

nexttile
pspectrum(v_data_filt, data_fs,'persistence',...
    'FrequencyLimits',[f_min f_lim],'TimeResolution',t_res)
title('Raw Data with sgolay Filter')
set(gca,'FontSize',font_size)

nexttile
pspectrum(v_right_interp,sim_fs,'persistence',...
    'FrequencyLimits',[f_min f_lim],'TimeResolution',t_res)
title('Simulation Data')
set(gca,'FontSize',font_size)

%% Spectrograms

f_lim = 300000;
t_res = 0.1;
% leakage = 0.85;
close all
tiledlayout(3,1);
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.58 0.9],'MenuBar','none');
nexttile
pspectrum(v_data, data_fs,'spectrogram',...
    'FrequencyLimits',[0 f_lim],'TimeResolution',t_res,'Leakage',0.5)
title('Unfiltered Raw Data')
set(gca,'FontSize',font_size)

nexttile
pspectrum(v_data_filt, data_fs,'spectrogram',...
    'FrequencyLimits',[0 f_lim],'TimeResolution',t_res,'Leakage',0.5)
title('Raw Data with sgolay Filter')
set(gca,'FontSize',font_size)

nexttile
pspectrum(v_right_interp,sim_fs,'spectrogram',...
    'FrequencyLimits',[0 f_lim],'TimeResolution',t_res,'Leakage',0.85)
title('Simulation Data')
set(gca,'FontSize',font_size)


%% Animation of Strain/Stress vs. Velocity
r = linspace(0,1,N+2)';
g = zeros(N+2,1);
b = linspace(1,0,N+2)';
rgb = [r,g,b];

tspan = 1:20000;
dx = N/L_rod;
k = (1:N)*2;
t_sim_norm = t_sim*c/L_rod;

% Redefine vars over smaller time span
U = S.y(k-1,tspan);
DU = diff(U);
DT = diff(t_sim(tspan));
DUDT = diff(U,1,2)./DT';
DUDX = DU./dx;
DUDX(N,:) = zeros(1,length(tspan));

%% Check Alignment of Waveforms
period_shift = 1; % starting point of animation (~ 440 periods = 100 msec)
offset = 0; % small shift due to drifting of waveform

t_range = [0 4]; % Specify width of plots in the time dimension
period_shift_msec = (t_range+period_shift+offset)*tau*1000;
t_strain_offset = t_sim_norm - period_shift - offset - L_rod; 

close all
tiledlayout(2,1);
set(gcf,'color','w');
nexttile
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.9],'MenuBar','none');
plot(t_sim(tspan(1):tspan(end-1))*1000 - L_rod, DUDT(1,:), 'LineWidth', 3, 'Color', rgb(1,:));
grid on
title('Velocity vs. Time')
xlabel('Time (msec)')
ylabel('Displacement (m)')
xlim(period_shift_msec/2);
set(gca,'FontSize',20)
% ylim([-0.02 0.1]);

nexttile
plot(t_strain_offset(tspan), DUDX(1,:), 'LineWidth', 3, 'Color', 'r');
xlim(t_range)
% ylim([-0.020 0.1]);
set(gca,'FontSize',20)
xlabel('Position Along Rod');
ylabel('Strain')
title('Strain vs. Time Normalized')
grid on

%% Animate Velocity/Strain (displayed time is off by 2x)
n_periods = 15;
period_shift = 0;%10000; %883 0.13
offset = 0; %0.5;
t_range = [0 1]; % Specify width of plots in the time dimension

% F(N*n_periods) = struct('cdata',[],'colormap',[]); 
% v = VideoWriter('StrainWaveAlongRod100Mfrom500msec', 'MPEG-4');
% v.FrameRate = round((3*N)/10);
% open(v)

t_vel_offset = t_sim_norm - period_shift - offset - L_rod; 
t_strain_offset = t_sim_norm - period_shift - offset - L_rod;

tau_i = 0;
tau_i = tau_i + period_shift/2; % for displaying time  during animation
dim = [0.46 0.42 0.2 0.1];
fig_size = [0.44 0.05 0.56 0.9];

close all
figure(1)
set(gcf,'color','w');
fig = gcf;
subplot(2,1,1)
set(gcf,'units','normalized','position',fig_size,'MenuBar','none');
p1 = plot(t_vel_offset(tspan(1):tspan(end-1)), DUDT(1,:), 'LineWidth', 3, 'Color', 'b');
set(p1, 'xdata', t_vel_offset(tspan(1):tspan(end-1)), 'ydata', DUDT(1,:));
ax1 = gca;
ax1.XLim = t_range;
ax1.YLim = [-0.001 0.01];%[0.0025 0.0032];
set(gca,'FontSize',20)
ax1.XLabel.String = 'Position Along Rod';
ax1.YLabel.String = 'Velocity [m/s]';
title('Velocity vs. Time')
grid on

subplot(2,1,2)
p2 = plot(t_strain_offset(tspan), DUDX(1,:), 'LineWidth', 3, 'Color', 'b');
set(p2, 'xdata', t_strain_offset(tspan), 'ydata', DUDX(1,:));
ax2 = gca;
ax2.XLim = t_range;
ax2.YLim = [-0.35*10^-10 0.35*10^-10];%[-0.02*10^-10 0.02*10^-10];
set(gca,'FontSize',20)
ax2.XLabel.String = 'Position Along Rod';
ax2.YLabel.String = 'Strain';
str = ['Time: ',num2str(round(tau*tau_i*1000,1)),' msec'];
text = annotation(fig,'textbox',dim,'String',str,'FitBoxToText','on');
text.FontSize = 28;
title('Strain vs. Time')
grid on


for i = 1:n_periods*N
    if mod(i,N) == 0 && i ~= 0
        tau_i = tau_i + 0.5;
        text.String = ['Time:  ', num2str(round(tau*tau_i*1000,1)),'  msec'];
        t_vel_offset = t_vel_offset - 1; 
        t_strain_offset = t_strain_offset - 1; 
    end
    n = mod(i,N)+1;
    set(p1, 'xdata', t_vel_offset(tspan(1):tspan(end-1)), 'ydata', DUDT(n,:));
    set(p2, 'xdata', t_strain_offset(tspan), 'ydata', DUDX(n,:));  
    % F(i) = getframe(gcf);
    % writeVideo(v, F(i))
    pause(0.01)
end
% close(v)

%% 3-D Experiments: Check Alignment of Waveforms
period_shift = 1; % starting point of animation (~ 440 periods = 100 msec)
offset = 0; % small shift due to drifting of waveform

t_range = [0 4]; % Specify width of plots in the time dimension
period_shift_msec = (t_range+period_shift+offset)*tau*1000;
fig_size = [0.44 0.05 0.56 0.9];

dt = tspan(1):tspan(end-1);

close all
figure(1)
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot3(t_sim(dt)*1000 - L_rod, DUDT(1,:),DUDX(1,1:end-1),'LineWidth', 3, 'Color', rgb(1,:));
grid on
title('Velocity vs. Time')
xlabel('Time (msec)')
ylabel('Displacement (m)')
xlim(period_shift_msec/2);
set(gca,'FontSize',20)
view(-10,15)
% ylim([-0.02 0.1]);

%% 3-D Experiments Animation
n_periods = 15;
period_shift = 0;%10000;
% period_shift = 883; %883 0.13
offset = 0; %0.5;
t_range = [0 2]; % Specify width of plots in the time dimension

t_vel_offset = t_sim_norm - period_shift - offset - L_rod; 
t_strain_offset = t_sim_norm - period_shift - offset - L_rod;

tau_i = 0;
tau_i = tau_i + period_shift/2; % for displaying time  during animation
dim = [0.46 0.42 0.2 0.1];
fig_size = [0.44 0.2 0.56 0.7];
fig_view = [25,15];
dt = tspan(1):tspan(end-1);

close all
figure(1)
set(gcf,'color','w');
fig = gcf;
set(gcf,'units','normalized','position',fig_size,'MenuBar','none');
p1 = plot3(t_vel_offset(dt), DUDT(1,:),DUDX(1,1:end-1)*10^11, 'LineWidth', 6, 'Color', 'b');
set(p1, 'xdata', t_vel_offset(dt), 'ydata', DUDT(1,:),'zdata',DUDX(1,1:end-1)*10^11);
ax1 = gca;
ax1.XLim = t_range;
ax1.YLim = [-0.001 0.01];%[0.0025 0.0032];
ax1.ZLim = [-3.5 3.5];
set(gca,'FontSize',20)
ax1.XLabel.String = 'Position Along Rod';
ax1.YLabel.String = 'Velocity [m/s]';
ax1.ZLabel.String = 'Strain';
title('Velocity vs. Time')
view(fig_view)
grid on

for i = 1:n_periods*N
    if mod(i,N) == 0 && i ~= 0
        %tau_i = tau_i + 0.5;
        %text.String = ['Time:  ', num2str(round(tau*tau_i*1000,1)),'  msec'];
        t_vel_offset = t_vel_offset - 1; 
    end
    n = mod(i,N)+1;
    set(p1, 'xdata', t_vel_offset(dt), 'ydata', DUDT(n,:),'zdata',DUDX(n,1:end-1)*10^11);
    view(fig_view)
    pause(0.01)
end












