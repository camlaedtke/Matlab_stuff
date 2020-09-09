%% CHAOTIC CIRCUIT
close all
clear all

% IMPORTANT PARAMETERS
Fs = 60; % Sampling frequency 
Ns = 18000; % Number of samples
t_end = Ns/Fs;
noisepow = -100; % dBW
time = t_end; % seconds
V_0 = 0.25; % V
R_0 = 157000; % R_0 = 157; % Ohms
R = 47000; % R = 47;
C = 1*10^-6; % Farads
R_1 = 12000; % R_1 = 150; 
R_2 = 82000; % R_2 = 889; 
R_nu = 82000;% Chaos between ~65k and ~110k

% Audible range
% Fs = 20000; % Sampling frequency 
% Ns = 20000; % Number of samples
% t_end = Ns/Fs;
% noisepow = -70; % dBW
% time = t_end; % seconds
% V_0 = 0.25; % V
% R_0 = 157; % Ohms
% R = 47;
% C = 1*10^-6; % Farads
% R_1 = 150; 
% R_2 = 889; 
% R_nu = 102;% Chaos between ~65 and ~110

% gamma,beta,omega aren't necessary, but they clean things up
gamma = R_2/R_1; % controls type of chaotic behavior
beta = R/R_nu; % controls chaos
omega = R/R_0; % controls amplitude of signal
t_f = time/(R*C); % dimensionless time

fprintf('\nSimulating over %.0f seconds \n',time)
fprintf('Characteristic Frequency f = %.2f Hz\n',1/(2*pi*R*C))

% Initial conditions
Y_0 = [V_0,0,0]; 
tspan = [0, t_f];

% Function handle
M = @(t,Y)[Y(2);Y(3);...
          (-omega - beta.*Y(3)-gamma.*Y(1) - Y(2))*(Y(1)<0) + ...
          (-omega - beta.*Y(3)-Y(2))*(Y(1)>0)];
options = odeset('Stats','on','RelTol',1e-13,'AbsTol',1e-13);
S = ode45(M,tspan,Y_0, options);

t = S.x';
tstep = round(length(t)/time);


fprintf('Time: %.0f steps per second\n',tstep)

plotstart = 1;% plotstart = round(1*tstep);
t = t(plotstart:end);% Cut out first few seconds of wacky behavior
t = t.*(R*C); % convert dimensionless time back into regular time
y = S.y';
V = y(plotstart:end,1);
dV = y(plotstart:end,2);
ddV = y(plotstart:end,3);

% voltage vs time
figure(1)
set(gcf,'units','normalized','position',[0.7 0.02 0.3 0.3],'MenuBar','none');
scatter(t,V,3,'red')
hold on
plot(t,V,'b')
plot(t,dV,'b','LineStyle',':','LineWidth',1.5)
plot(t,ddV,'m','LineStyle',':','LineWidth',1.5)


findpeaks(V,t) % plot peaks
[pks,locs] = findpeaks(V,t); % get peaks
title('V(t) Time Domain')
xlabel('Time (seconds)')
ylabel('V')
xlim([locs(10) locs(15)])
set(gca,'FontSize',15)
grid on
hold off

% 2-D phase space
figure(2)
set(gcf,'units','normalized','position',[0.4 0.64 0.3 0.3],'MenuBar','none');
plot(V,dV,'b')
title('2-D Phase Space')
xlabel('$V$','interpreter','latex')
ylabel('$\partial V / \partial t$','interpreter','latex')
set(gca,'FontSize',15)
grid on

% 3-D phase space
figure(3)
set(gcf,'units','normalized','position',[0.7 0.64 0.3 0.3],'MenuBar','none');
plot3(V,dV,ddV,'b','LineWidth',0.3)
title('3-D Phase Space','interpreter','latex')
xlabel('$V$','interpreter','latex')
ylabel('$\partial V / \partial t$','interpreter','latex')
zlabel('$\partial^{2}V/\partial t^2$','interpreter','latex')
set(gca,'FontSize',15)
grid on
view(30,30)


% FREQUENCY DOMAIN
f_max = round(Fs/(2)); % Max frequency (for plots)
res = Fs/Ns; % Frequency resolution (for power spectrum plot)

Y = fft(V); % Fourier Transform
P = abs(Y*2/Ns); % Realize and Normalize
f = Fs*(0:Ns-1)/Ns; % Frequency list in Hz (bin widths)

fprintf('Number of Samples: %.0f \n',Ns)
fprintf('Sampling Frequency: %.0f Hz \n',Fs)
fprintf('Frequency Resolution: %.3f Hz \n',res*10^3)
res = res+0.01; % it rounds down sometimes, which causes an error

figure(4)
P(1) = 0;
set(gcf,'units','normalized','position',[0.4 0.02 0.3 0.3],'MenuBar','none');
scatter(f(1:(Ns/2)),P(1:(Ns/2)),3,'red')
hold on
plot(f(1:Ns/2),P(1:Ns/2),'b')
title('Single-Sided Amplitude Spectrum of V(t)','interpreter','latex')
xlabel('Frequency (Hz)','interpreter','latex') 
ylabel('$log \left( |P(f)| \right)$','interpreter','latex')
set(gca,'yscale','log')
set(gca,'FontSize',15)
grid on
hold off

V = V - mean(V); % Remove DC Offset (good for power spectra)

% Power spectrum
xTable = timetable(seconds(t),V); 
figure(5)
set(gcf,'units','normalized','position',[0.4 0.33 0.3 0.3],'MenuBar','none');
[pxx,f] = pspectrum(xTable,'FrequencyLimits',[0 f_max],...
    'FrequencyResolution',10*res);
plot(f,pow2db(pxx),'b','LineWidth',1)
title('Power Spectrum','interpreter','latex')
xlabel('Frequency (Hz)','interpreter','latex')
ylabel('Power (dBW)','interpreter','latex')
ylim([noisepow 0])
set(gca,'FontSize',15)
grid on

% Spectrogram
figure(6)
set(gcf,'units','normalized','position',[0.7 0.33 0.3 0.3],'MenuBar','none');
[sp,fp,tp] = pspectrum(xTable,'spectrogram','FrequencyLimits',[0 f_max],...
    'FrequencyResolution',10*res);
sp = pow2db(sp);
tp = seconds(tp);
mesh(tp,fp,sp)
view(110,25)
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlim([noisepow 0])
zlabel('Power (dBW)')
set(gca,'FontSize',15)






