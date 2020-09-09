
clear all

% IMPORTANT PARAMETERS
Fs = 30; % Sampling frequency 
Ns = 9000; % Number of samples
t_end = Ns/Fs;
noisepow = -100; % dBW
time = t_end; % seconds
V_0 = 0.25; % V
R_0 = 157000; % R_0 = 157; % Ohms
R = 47000; % R = 47;
C = 1*10^-6; % Farads
R_1 = 15000; % R_1 = 150; 
R_2 = 89000; % R_2 = 889; 
R_nu = 55000;% Chaos between ~65k and ~110k

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

figure(1)
set(gcf,'units','normalized','position',[0.55 0.1 0.4 0.8],'MenuBar','none');
scatter(t,V,3,'red')
hold on
plot(t,V,'b')
findpeaks(V,t) % plot peaks
[pks,locs] = findpeaks(V,t); % get peaks
title('Voltage vs Time','interpreter','latex')
xlabel('Time','interpreter','latex')
ylabel('Voltage','interpreter','latex')
xlim([locs(150) locs(165)])
ylim([-0.3 3.5])
set(gca,'FontSize',32)
grid on
hold off

% figure(2)
% P(1) = 0;
% set(gcf,'units','normalized','position',[0.4 0.33 0.6 0.3],'MenuBar','none');
% scatter(f(1:(Ns/2)),P(1:(Ns/2)),3,'red')
% hold on
% plot(f(1:Ns/2),P(1:Ns/2),'b')
% title('Single-Sided Amplitude Spectrum of V(t)')
% xlabel('Frequency (Hz)') 
% ylabel('|P(f)|')
% set(gca,'FontSize',15)
% grid on
% hold off



