% FULL SCRIPT

clear all
close all
DataFileName ={'tek0039All.mat';'tek0033All.mat'};
RodLength =[0.551;0.399];
RodDiameter = [0.01271;0.01271];
MyData = table(DataFileName, RodLength,RodDiameter);
% longer rod data
load(MyData.DataFileName{1})

%%
close all
% Leading edge of trigger signal: t impact = 0.0048703 seconds
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.35 0.3 0.65 0.44]);
plot(t,X,'LineWidth',1)
hold on
plot(t,Y,'LineWidth',1)
plot(t,V,'LineWidth',4)
grid on
title('Light Intensity at Detectors X and Y and Trigger Signal vs. Time')
xlabel('Time [seconds]','Interpreter','latex')
ylabel('Voltage [V]','Interpreter','latex')
legend('VX', 'VY', 'Trigger Signal')
set(gca,'FontSize',28)
xlim([0.00486 0.00502])
t_end = 0.004999;
t_start = 0.00487;
t_contact = t_end - t_start;
fprintf('\nContact time: t_c = %.4e\n',t_contact)

%% Converting Intensity Data to Displacements at the End of the Rod
close all
Xmedian  = median(X);
Ymedian  = median(Y);

%calculate dphi and unwrap it
phiraw =atan2(Y-Ymedian, X-Xmedian);
phirawunw = unwrap(phiraw);

%convert to displacement
lam = 632.8E-9;  %HeNe laser wavelength
du = -phirawunw*lam/(pi*4); %calculate the displacement
% dudt = diff(du)./diff(t);

% save('tek0039RawDisplacementAndVelocity.mat','du','dudt','t')
%% Mess around with the raw data before its averaged
close all

font_size = 18;
view_shift = [160,10];

dt = 1:1:length(t)-2;
du_filt = sgolayfilt(du,5,91);
dudt = diff(du_filt)./diff(t);
d2udt2 = diff(dudt)./diff(t(1:end-1));

figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56,0.44],'MenuBar','none');
plot(t(1:length(t)-1)*10^3, 10^6*(diff(du)./diff(t)),'LineWidth', 3, 'Color', 'b')
hold on
plot(t(dt)*10^3, 10^6*dudt(dt),'LineWidth', 3, 'Color', 'g')
hold off
grid on
title(' Velocity Data')
xlabel('Time [msec]','Interpreter','latex')
xlim([5.9 6])
ylabel('Velocity [um/sec]','Interpreter','latex')
set(gca,'FontSize',font_size)

[velFT,f,Fs] = MyFFT(dudt, t); 
figure(3)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
plot(f*10^-3, velFT,'LineWidth', 1, 'Color', 'k')
grid on
title('Fourier Transform of Velocity Data diff(solayfilt(du,5,91))./diff(t))')
xlabel('Frequency [kHz]','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')
xlim([0 1000])
ylim([0.00002*10^-6 2*10^-2])
set(gca,'FontSize',font_size,'YScale', 'log')


%% NOISE AND DATA REDUCTION
% noise amplitude on the data is about 1 nm, so we average over 20 dps
blocksize = 20;   
newrows = floor(length(t)/blocksize);
t_ave_full = mean(reshape(t(1:newrows*blocksize), blocksize, newrows))'; 
du_ave_full = mean(reshape(du(1:newrows*blocksize), blocksize, newrows))'; 
V_ave = mean(reshape(V(1:newrows*blocksize), blocksize, newrows))'; 
S_ave = mean(reshape(S(1:newrows*blocksize), blocksize, newrows))'; 

vel_ave = diff(du_ave_full)./diff(t_ave_full);
vel_ave(end+1) = 0;

 % only use a subset of the data, i.e., the first 10 msec
dt = 1: 0.01*length(t_ave_full);

close all
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.34 0.44]);
plot(t*10^3, du*10^9,'LineWidth', 2, 'Color', 'b','DisplayName','Raw Data')
hold on
plot(t_ave_full*10^3, du_ave_full*10^9,'LineWidth', 3, 'Color', 'r','DisplayName','Averaged')
grid on
title('Displacement vs. Time')
xlabel('Time [msec]','Interpreter','latex')
xlim([5.15 5.155])
ylabel('Displacement [nm]','Interpreter','latex')
legend
set(gca,'FontSize',22)

figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.7 0.5 0.34 0.44]);
plot(t_ave_full*10^3, du_ave_full*10^3,'LineWidth', 3, 'Color', 'r','DisplayName','Averaged')
grid on
title('Displacement vs. Time')
xlabel('Time [msec]','Interpreter','latex')
xlim([4.5 6])
ylabel('Displacement [mm]','Interpreter','latex')
legend
set(gca,'FontSize',22)

length(dt)
% save('tek0039Averaged.mat','t_ave_full','du_ave_full','vel_ave','V_ave','S_ave','dt')

%% GET CM VELOCITY , PARTICLE VELOCITY, PHASE VELOCITY, AND FREQUENCY
clear all
close all
load('tek0039Averaged.mat') % all relevant variables after averaging
disp('-------------- Preliminary Analysis ---------------')

% Toss out all data before the impact – it’s not relevant.
t_impact = 0.0048703; % rough estimate based on trigger signal
idx_impact = find(abs(t_ave_full(dt)-t_impact) < 0.000001); %impact idx
idx_impact = idx_impact (1); % select the first one (if there are multiple)
t_ave = t_ave_full - t_impact; % shift so zero is at impact time
dt = idx_impact:(idx_impact+4999); % make dt start at impact time 

%% Filter first 10 msec of displacement data

font_size = 22;
% t_span = t_ave(dt(end)) - t_ave(dt(1));
% tau = 2.2619e-4;
% n_periods = t_span/tau;

fs = 1/(t_ave(2)-t_ave(1));

du_ave_lowpass = lowpass(du_ave_full(dt),4423,500000);
% du_ave_lowpass = sgolayfilt(du_ave_full(dt),1,45);
% du_ave_filter = filter(du_ave_full(dt))

% Least-squares-fit the first 10 msec to get com velocity 
[fitobj, gof, outp] = fit(t_ave(dt), du_ave_lowpass,'poly1');
modelparams.v_cm_fit = fitobj;

close all
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44]);
p1 = plot(t_ave(dt), fitobj(t_ave(dt)*10^3),'LineWidth', 3, 'Color', 'r');
hold on
% plot(t_ave(dt), du_ave_lowpass*10^3,'LineWidth', 5, 'Color', 'g','DisplayName','Lowpass')
p2 = plot(t_ave(dt), du_ave_full(dt)*10^3,'LineWidth', 3, 'Color', 'k');
hold off
grid on
title('LSQ Fit of First 10 msec of Displacement Data')
xlabel('Time (seconds)')
xlim([0 0.01])
ylim([0 3e-1])
ylabel('Displacement (mm)')
legend('LSQ Fit','Data','Location','northwest')
set(gca,'FontSize',font_size)


% figure(2)
% set(gcf,'color','w');
% set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
% p1 = plot(t_ave(dt), fitobj(t_ave(dt)*10^3),'LineWidth', 3, 'Color', 'r');
% hold on
% % plot(t_ave(dt), du_ave_lowpass,'LineWidth', 5, 'Color', 'g','DisplayName','Lowpass')
% p2 = plot(t_ave(dt), du_ave_full(dt)*10^3,'LineWidth', 3, 'Color', 'k');
% hold off
% grid on
% title('LSQ Fit of First 10 msec of Displacement Data')
% xlabel('Time (seconds)')
% xlim([0.004 0.006])
% % ylim([0 3e-1])
% ylabel('Displacement (mm)')
% legend('LSQ Fit','Data','Location','northwest')
% set(gca,'FontSize',font_size)

figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
plot(t_ave(dt), vel_ave(dt),'LineWidth', 3, 'Color', 'r');
hold off
grid on
title('Velocity of Rod at Free End')
xlabel('Time (seconds)')
xlim([0 0.001])
% ylim([0 3e-1])
ylabel('Velocity (m/s)')
set(gca,'FontSize',font_size)

% Obtain Parameter Values
v_cm = fitobj.p1; % Slope of displacement vs time
v_cm_err = gof.rmse;
du_0 = du_ave_full(1); % initial displacement
du_ave = du_ave_full - du_0; % shift so zero is at initial displacement
y_0 = fitobj.p2; % y-intercept of v cm fit
fprintf(['Center of Mass Velocity \n     v_cm = (%.3f +/- %.3f)*10^-3  m/s\n '...
        '    redchi = %.4f \n'],v_cm*10^3, gof.rmse*10^3, gof.adjrsquare)

% Get characteristic frequency
[velFT,f,Fs] = MyFFT(vel_ave, t_ave); 
[~, f0] = findpeaks(velFT, f,'MinPeakHeight',0.008);
f0 = f0(2);

figure(3)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.03 0.56 0.44],'MenuBar','none');
plot(f, velFT,'LineWidth', 3, 'Color', 'b')
hold on
findpeaks(velFT, f,'MinPeakHeight',0.008);
hold off
grid on
title('Velocity vs Time Fourier Transform')
xlabel('Frequency (Hz)')
ylabel('Velocity (m/s)')
xlim([0 20000])
set(gca,'FontSize',font_size);%,'YScale', 'log')

%  Uncertainty is fs/N
f0_err = Fs/length(vel_ave); % Hz
w0 = 2*pi*f0; % rad/sec
w0_err = 2*pi*f0_err;
fprintf(['Characteristic Frequency \n'...
 '     w0 = %.0f +/- %.0f rad/sec   \n'...
 '     f0 = %.0f +/- %.0f Hz \n'],w0,w0_err,f0,f0_err);
tau = 1/f0;
tau_err = sqrt(f0_err*(1/f0^2)^2);
fprintf('     tau = %.4e +/- %.2e seconds \n',tau,tau_err)


% Get particle velocity
v_particle = 2*v_cm; % particle velocity (over dt)
v_particle_err = sqrt(v_cm_err^2*(2)^2);
fprintf(['Particle Velocity \n'...
 '     v_particle = (%.3f +/- %.3f)*10^-3 m/s\n'],v_particle*10^3, v_particle_err*10^3);

% Dimensions of rod...
L_rod = 0.551; % meters
L_rod_err = 0.001;
D_rod = 0.01271; % meters
D_rod_err = 0.00001; % meters
A_rod = pi*(D_rod/2)^2; % meters^2
V_rod = A_rod*L_rod;
V_rod_err = sqrt(L_rod_err^2*(pi*D_rod^2/4)^2+D_rod_err^2*(pi*D_rod*L_rod/2)^2);

% Dimensions of ball ...
D_ball = 0.02404; % meters
D_ball_err = 0.0001; % meters
V_ball = (pi/6)*D_ball^3;
V_ball_err = sqrt(D_ball_err^2*(pi*D_ball^2/2)^2);

% For 304?? stainless steel ...
rho = 8000; % kg/m^2 
rho_err = 100; 

% Get phase velocity (c0)
v_phase = (2*L_rod)/tau;
v_phase_err = sqrt(L_rod_err^2*(2/tau)^2+tau_err^2*(2*L_rod/tau^2)^2);
fprintf(['Phase Velocity \n'...
 '     v_phase = %.0f +/- %.0f m/s\n'],v_phase,v_phase_err);

% Get Young's Modulus
Y = rho*v_phase^2; % Pa
Y_err = sqrt(rho_err^2*(v_phase)^2 + v_phase_err^2*(rho)^2);
fprintf('Youngs Modulus \n     Y = %.4f +/- %.4f GPa\n',Y*10^-9, Y_err*10^-9);

% Get mass of rod
m_rod = rho*V_rod; % kg
m_ball = rho*V_ball; %kg

fprintf('Dimensions of Ball and Rod \n')
fprintf('     V rod = (%.3f +/- %.3f)*10^-6 m^3\n',V_rod*10^6,V_rod_err*10^6)
fprintf('     V ball = (%.3f +/- %.3f)*10^-6 m^3\n',V_ball*10^6,V_ball_err*10^6)
fprintf('     m_rod = %.4f kg\n',m_rod);
fprintf('     m_ball = %.4f kg\n\n',m_ball);

save('tek0039Params.mat','dt','t_ave','du_ave','vel_ave',...
    "v_phase","v_particle","v_cm","w0","L_rod",'D_rod','D_ball','rho','A_rod','V_rod','Y')
save('tek0039ImportantParams.mat','v_phase','v_particle',...
    'm_rod','m_ball','tau','w0','f0','Y','L_rod','D_ball','A_rod')

%% Remove Motion of Pendulum

t_post_impact = t_ave(idx_impact:length(t_ave)); % toss out pre-impact times
du_post_impact = du_ave(idx_impact:length(t_ave)); % toss out pre-impact du
    
% STEP 2:  to fit the raw displacement data to a 3rd degree polynomial 
% over 20001 datapoints to obtain the (lowpass filtered) background motion
order = 3;
framelen = 20001;
du_sgolay = sgolayfilt(du_post_impact, order, framelen);

font_size = 20;
close all
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44],'MenuBar','none');
plot(t_post_impact , du_sgolay, 'LineWidth', 3, 'Color', 'g')
hold on
plot(t_post_impact , du_post_impact,'--', 'LineWidth', 1, 'Color', 'k')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([t_post_impact(105000) t_post_impact(115000)]) 
grid on
set(gca,'FontSize',16)

% STEP 3: Subtract the displacement data (point-by-point) from the 
% corresponding filtered background data from step 2.
du_pendulum_subtracted = du_post_impact - du_sgolay;

% Figures represents displacement data we would see if we were moving
% with the center of mass of the rod.
figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.28 0.44],'MenuBar','none');
plot(t_post_impact, du_pendulum_subtracted, 'LineWidth', 0.1, 'Color', 'b')
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
grid on
set(gca,'FontSize',16)

% STEP 4: add the center of mass displacement calculated from the
% linear fit from step 1 back to the data,
du_cm = v_cm*t_post_impact + y_0;
du_no_gravity = du_pendulum_subtracted + du_cm;
du_no_gravity = du_no_gravity - du_no_gravity(1);

% Figures represents displacement data we would see if we were moving
% with the center of mass of the rod.
figure(3)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.7 0.05 0.28 0.44],'MenuBar','none');
plot(t_post_impact, du_no_gravity, 'LineWidth', 2, 'Color', 'b')
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
grid on
set(gca,'FontSize',16)

v_no_gravity = diff(du_no_gravity)./diff(t_post_impact);
v_no_gravity(end+1) = 0;


%%
save('tek0039PendulumRemoved.mat','v_cm','t_post_impact','du_post_impact',...
    'du_pendulum_subtracted','du_no_gravity','v_no_gravity')
