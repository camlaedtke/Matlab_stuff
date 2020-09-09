% Ball Rod Priliminary % Loading the Data
clear all;
close all
DataFileName ={'tek0039All.mat';'tek0033All.mat'};
RodLength =[0.551;0.399];
RodDiameter = [0.01271;0.01271];
MyData = table(DataFileName, RodLength,RodDiameter)
% longer rod data
load(MyData.DataFileName{1})

%% Raw Data Plotting
close all
disp('Figure 3')

figure(1)
set(gcf,'color','w','MenuBar','none');
set(gcf,'units','normalized','position',[0.42 0.5 0.58 0.44]);
plot(t, X, t, Y, t, V)
grid on
title('Light Intensity at Detectors X and Y and Trigger Signal vs. Time')
xlabel('Time (seconds)')
ylabel('Voltage')
legend('gx IX', 'gy IY', 'Trigger Signal')
set(gca,'FontSize',22)
% xlim([0.09 0.095])

% we can observe from Figure 3 that the impact (for this data set) 
% occurred at around t = 5 msec.

disp('Figure 4')
figure(2)
set(gcf,'color','w','MenuBar','none');
set(gcf,'units','normalized','position',[0.42 0.05 0.58 0.44]);
plot(t, X, t, Y, t, V)
grid on
title('Light Intensity at Detectors X and Y and Trigger Signal vs. Time')
xlabel('Time (seconds)')
ylabel('Voltage')
legend('gx IX', 'gy IY', 'Trigger Signal')
set(gca,'FontSize',22)
xlim([0.0049 0.0054])


%% Converting Intensity Data to Displacements at the End of the Rod
close all
% Plotting the voltages from the two detectors against each other
disp('Figure 5')
figure(3)
set(gcf,'color','w','MenuBar','none');
set(gcf,'units','normalized','position',[0.5 0.2 0.55 0.75]);
plot(X, Y,'LineWidth', 0.3)
grid on
title('Raw Data: Intensity at Detector X vs. Detector Y')
xlabel('Voltage at Detector X')
ylabel('Voltage at Detector Y')
set(gca,'FontSize',24)

% find the "center" of circular data, Xmedian and Ymedian and then subtract
% it.
Xmedian  = median(X);
Ymedian  = median(Y);

%calculate dphi and unwrap it
phiraw =atan2(Y-Ymedian, X-Xmedian);
phirawunw = unwrap(phiraw);

%convert to displacement
lam = 632.8E-9;  %HeNe laser wavelength
du = -phirawunw*lam/(pi*4); %calculate the displacement

%% Preliminary Look at the Displacement Data
close all

% The displacement of the rod?s end, opposite its impact point, 
% over the entire 1 second time interval over which data was acquired.

disp('Figure 6')
figure(4)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.3 0.44]);
plot(t, du*1000, 'LineWidth', 2.5)
grid on
title('Displacement vs. Time','FontWeight','bold','Color','b' )
ylabel('Distance (mm)','FontWeight','bold')
xlabel('Time (sec)','FontWeight','bold')
set(gca,'FontSize',22)

% We can clearly observe the compressional waves when we zoom in 
% to the first 10 msec of the data acquisition as shown in Figure 7. 
% (Compare this to Figure 3, where we show the intensity and 
% trigger data over the same time interval.)

disp('Figure 7')

figure(5)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
plot(t, du, 'LineWidth', 2)
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0.0 0.01])
set(gca,'FontSize',22)

% As we further zoom in, see Figure 8, and only study the first 3 
% compression waves (by limiting the time axis to a 600 usec interval) 
% we can clearly observe on the 2nd plateau, 5.3 msec < t < 5.4 msec, 
% that the end of the rod reversed its direction.

disp('Figure 8')

figure(6)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.28 0.44],'MenuBar','none');
plot(t, du, 'LineWidth', 2)
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0.0049 0.0055]) 
set(gca,'FontSize',22)

vel = diff(du)./diff(t);

disp('Figure 8')
figure(7)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.7 0.05 0.28 0.44],'MenuBar','none');
plot(t(1:end-1),vel, 'LineWidth', 2)
grid on
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0.0049 0.0055]) 
set(gca,'FontSize',22)


%% NOISE AND DATA REDUCTION
close all

%reduce the data by averaging over blocks of data points
blocksize = 20;    %datapoints to average over (17?)

newrows = floor(length(t)/blocksize); % calculate the (integer) number of rows of the resulting matrix
t_ave = mean(reshape(t(1:newrows*blocksize), blocksize, newrows))'; %average the time data
du_ave = mean(reshape(du(1:newrows*blocksize), blocksize, newrows))'; %average the displacement data
V_ave = mean(reshape(V(1:newrows*blocksize), blocksize, newrows))'; %average the trigger data
S_ave = mean(reshape(S(1:newrows*blocksize), blocksize, newrows))'; %average the audio data

disp('Figure 9')
figure(7)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44],'MenuBar','none');
plot(t, 1E6*du,'Color','k', 'LineWidth', 2)
hold on
plot(t_ave, 1E6*du_ave,'g', 'LineWidth', 1)
% plot(t, 1E6*sgolayfilt(du,3,21),'g','LineWidth', 2)
hold off
grid on
title('Displacement vs. Time')
xlabel('Time (seconds)','interpreter','latex')
ylabel('Displacement (um)','interpreter','latex')
xlim([0.0051 0.0052])
% legend('Original Data','Averaged Data','sgolayfilt')
set(gca,'FontSize',22)

vel_sg = diff(sgolayfilt(du,7,51))./diff(t);
vel = diff(du)./diff(t);

figure(8)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
plot(t(1:end-1),vel,':','Color','k', 'LineWidth', 1)
hold on
plot(t(1:end-1), vel_sg,'Color','b','LineWidth', 3)
hold off
grid on
title('Velocity vs. Time')
xlabel('Time (seconds)','interpreter','latex')
ylabel('Displacement (um)','interpreter','latex')
xlim([0.005 0.00516])
% legend('Original Data','Averaged Data','sgolayfilt')
set(gca,'FontSize',22)

% [velFT,f] = MyFFT(vel_sg,t(1:end-1));  % calc.and plot its FFT

  % (Note: a constant (vertical)  offset was added to the audio data in 
% Figure 12 for illustration purposes.)


% figure(9)
% set(gcf,'units','normalized','position',[0.42 0.25 0.56 0.34],'MenuBar','none');
% plot(f, velFT, 'LineWidth', 1, 'Color', 'black')
% title('Filtered Velocity Data FFT')
% xlabel('Frequency')
% ylabel('Amplitude')
% grid on
% xlim([0, 500000])
% set(gca,'FontSize',22,'YScale','log')
% hold off
% 

%% AUDIO DATA

% close all
% Figure 10 displays a plot of the recorded audio signal 
% during the first 10 msec.  It covers the same time interval 
% as Figures 3 and 7.

disp('Figure 10')
dt = 1: 0.01*length(t_ave); % only use a subset of the data, i.e., the first 10 msec
figure(8)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44],'MenuBar','none');
plot(t_ave(dt), S_ave(dt),'LineWidth', 2,'Color', 'r')
title('Audio Signal vs. Time')
xlabel('Time (seconds)')
ylabel('Audio Signal (V)')
set(gca,'FontSize',22)
grid on

% Figure 11, covering the first 5 msec after the impact, 
% the displacement and the audio signals are overlaid.

disp('Figure 11')
figure(9)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44],'MenuBar','none');
plot(t_ave(dt), du_ave(dt), 'LineWidth', 2, 'Color', 'b')
hold on
plot(t_ave(dt), S_ave(dt)/200, 'LineWidth', 2, 'Color', 'r')
hold off
title('Reduced Audio Signal and Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement and Audio Signal (arbs)')
grid on
legend('Displacement','Audio Signal')
xlim([0.005, .01])
set(gca,'FontSize',22)


% FOURIER TRANSFORM OF THE AUDIO SIGNAL AND THE DISPlACEMENT VELOCITY

vel_ave = diff(du_ave);  %calculate the average displacement velocity
[velFT,f] = MyFFT(vel_ave(dt),t_ave(dt));  % calc.and plot its FFT

  % (Note: a constant (vertical)  offset was added to the audio data in 
% Figure 12 for illustration purposes.)


disp('Figure 12')
figure(10)
set(gcf,'units','normalized','position',[0.62 0.05 0.38 0.44],'MenuBar','none');
plot( f, velFT*100000, 'LineWidth', 2, 'Color', 'black')
hold on
[SFT,f] = MyFFT(S_ave(dt),t_ave(dt));  % Audio signal FFT
plot(f, SFT+0.0005, 'LineWidth', 2, 'Color', 'r')
hold off

title('Fourier Transform of the Audio Signal and Velocity')
xlabel('f (Hz)')
ylabel('Velocity and Audio Signal (arbs + arb offset)')
xlim([0, 20000])
ylim([0 5E-3])
grid on
legend('Velocity', 'Audio')
set(gca,'FontSize',22)

%%
close all
figure(10)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44],'MenuBar','none');
[SFT,f] = MyFFT(S_ave,t_ave);  % Audio signal FFT
plot(f, SFT, 'LineWidth', 1, 'Color', 'r')
hold off
title('Fourier Transform of the Audio Signal and Velocity')
xlabel('f (Hz)')
ylabel('Velocity and Audio Signal (arbs + arb offset)')
xlim([0, 20000])
% ylim([0 5E-3])
grid on
legend('Velocity', 'Audio')
set(gca,'FontSize',22,'YScale','log')


%%
close all

Fs = 5.0000e+05;
L = size(S_ave,1);
res = Fs/L;


% Spectrogram
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.3 0.56 0.44]);
[sp,fp,tp] = pspectrum(S_ave,Fs,'spectrogram','FrequencyResolution',300*res);
mesh(tp,fp,pow2db(sp))
view(100,20)
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([0 32000])
zlim([-130 -40])
zlabel('Power')
set(gca,'FontSize',15)



function [y, f] = MyFFT(x, t)
L = length(x);
y = fft(x);
% Compute the two-sided spectrum P2. 
% Then compute the single-sided spectrum P1 
% based on P2 and the even-valued signal length L.

P2 = abs(y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%Define the frequency domain f and calculate its range
P1(1)=0; %zero out DC offset
y = P1;
Fs = 1/(t(2)-t(1))
f = Fs*(0:(L/2))/L;
end


