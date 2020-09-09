%% PRELIMINARY ANALYSIS OF VOLTAGE DATA
close all
clear all
% IMPORTANT PARAMETERS
Fs = 30; % Hz
Ns = 9000; % samples
t_end = Ns/Fs; % end time
noisepow = -100; % dBW

% DATA INPUT
V = dlmread('VoltageData/LowR9ks30hz.csv',',',0,0); % V
t = linspace(0,t_end,Ns)'; %(seconds) vector of time values

% voltage vs time
figure(1)
set(gcf,'units','normalized','position',[0.7 0 0.3 0.3],'MenuBar','none');
scatter(t,V,25,'k')
hold on
findpeaks(V,t) % plot peaks
[pks,locs] = findpeaks(V,t); % get peaks
title('V(t) Time Domain')
xlabel('Time (seconds)')
ylabel('V')
xlim([locs(10) locs(20)])
set(gca,'FontSize',15)

% Spline interpolation
% query points for spline function (4 points between timesteps)
tq = linspace(0,t_end,Ns*4);
Vq = spline(t,V,tq);

% Spline plot
scatter(tq,Vq,5,'red')
plot(tq,Vq,'red','LineStyle','--','LineWidth',1.5)

% Derivative of spline plot
dVq = diff(Vq);
scatter(tq(1:end-1),dVq,5,'b')
plot(tq(1:end-1),dVq,'b','LineStyle',':','LineWidth',1.5)

% Second derivative
ddVq = diff(dVq);
scatter(tq(1:end-2),ddVq,5,'m')
plot(tq(1:end-2),ddVq,'m','LineStyle',':','LineWidth',1.5)

grid on
hold off

% 2-D phase space
figure(2)
set(gcf,'units','normalized','position',[0.4 0.6 0.3 0.3],'MenuBar','none');
plot(Vq(1:end-1),dVq,'b')
title('2-D Phase Space')
xlabel('$V$')
ylabel('$\partial V / \partial t$','interpreter','latex')
set(gca,'FontSize',15)
grid on

% 3-D phase space
figure(3)
set(gcf,'units','normalized','position',[0.7 0.6 0.3 0.3],'MenuBar','none');
plot3(Vq(1:end-2),dVq(1:end-1),ddVq,'b','LineWidth',0.3)
title('3-D Phase Space','interpreter','latex')
xlabel('$V$','interpreter','latex')
ylabel('$\partial V / \partial t$','interpreter','latex')
zlabel('$\partial^{2}V/\partial t^2$','interpreter','latex')
set(gca,'FontSize',25)
grid on
view(30,30)


% FREQUENCY DOMAIN

L = size(V,1); % Data length (Length of signal)
Y = fft(V); % Fourier Transform
P = abs(Y*2/L); % Realize and Normalize
f = Fs*(0:L-1)/L; % Make bin widths
f_max = round(Fs/(2)); % Max frequency (for plots)
res = Fs/L; % Frequency resolution (for power spectrum plot)

fprintf('Number of Samples: %.0f seconds \n',Ns)
fprintf('Sampling Frequency: %.0f seconds \n',Fs)
fprintf('Frequency Resolution: %.3f seconds \n',res*10^3)

res = res+0.01; % when you have 'res' as an input to a function
% it rounds down sometimes, which causes an error

figure(4)
P(1) = 0;
set(gcf,'units','normalized','position',[0.4 0 0.3 0.3],'MenuBar','none');
scatter(f(1:(L/2)),P(1:(L/2)),3,'red')
hold on
plot(f(1:(L/2)),P(1:(L/2)),'b')
set(gca,'yscale','log')
title('Single-Sided Amplitude Spectrum of V(t)','interpreter','latex')
xlabel('Frequency (Hz)','interpreter','latex') 
ylabel('$log \left( |P(f)| \right)$','interpreter','latex')
set(gca,'FontSize',25)
grid on
hold off

V = V - mean(V); % Remove DC Offset

% Power spectrum
figure(5)
set(gcf,'units','normalized','position',[0.4 0.3 0.3 0.3],'MenuBar','none');
xTable = timetable(seconds(t),V); 
[pxx,f] = pspectrum(xTable,'FrequencyResolution',10*res);
plot(f,pow2db(pxx),'b','LineWidth',1)
title('Power Spectrum','interpreter','latex')
xlabel('Frequency (Hz)','interpreter','latex')
ylabel('Power (dBW)','interpreter','latex')
ylim([noisepow 0])
set(gca,'FontSize',25)
grid on

% Spectrogram
figure(6)
set(gcf,'units','normalized','position',[0.7 0.3 0.3 0.3],'MenuBar','none');
[sp,fp,tp] = pspectrum(V,Fs,'spectrogram','FrequencyResolution',10*res);
mesh(tp,fp,pow2db(sp))
view(110,25)
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([0 f_max])
zlim([noisepow 0])
zlabel('Power')
set(gca,'FontSize',15)






