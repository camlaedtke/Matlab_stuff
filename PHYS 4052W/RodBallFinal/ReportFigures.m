clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')
load('tek0039Sim100MFreeEndDataFixedV0.mat')

t_data = t_post_impact; 
u_data = du_no_gravity; v_data = v_no_gravity;
c = v_phase;

%% Velocity and Displacement over first msec
SCALE_RESULTS = 10;

close all
font_size = 20;

% t = tiledlayout(1,2);
% % t.TileSpacing = 'compact';
% % t.Padding = 'compact';
% t.Title.String = 'Particle Velocity at Free End';
% t.Title.FontWeight = 'bold';
% t.Title.FontSize = 28;


fig_size = [0.42 0.5 0.56 0.3];
% nexttile
figure(1)
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3,'DisplayName','')
hold on
plot(t_sim*10^3, dudt_free,'Color','b','LineWidth',3)
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
title('Particle Velocity at Free End')
xlabel('Time [msec]')
ylabel('Velocity [m/s]')
xlim([0 1])
ylim([-0.015 0.1])
set(gca,'FontSize',font_size)


% nexttile
% plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
% hold on
% plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
% hold off
% grid on
% % title('Particle Velocity at Free End')
% xlabel('Time [msec]')
% % ylabel('Velocity [m/s]','interpreter','latex')
% xlim([100 100.52])
% ylim([-0.015 0.1])
% set(gca,'FontSize',font_size)


%% Particle Velocity 400 msec
close all
font_size = 20;

figure(1)
t = tiledlayout(1,4);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
t.Title.String = 'Particle Velocity at Free End';
t.Title.FontWeight = 'bold';
t.Title.FontSize = 28;


fig_size = [0.3 0.6 0.7 0.34];
nexttile
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
ylabel('Velocity [m/s]')
xlim([100 100.52])
ylim([-0.015 0.08])
set(gca,'FontSize',font_size)


nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([200 200.52])
ylim([-0.015 0.08])
set(gca,'FontSize',font_size)


nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([300 300.52])
ylim([-0.015 0.08])
set(gca,'FontSize',font_size)

nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3,'DisplayName','Model')
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3,'DisplayName','Data')
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([400 400.52])
ylim([-0.015 0.08])
set(gca,'FontSize',font_size)
legend

% Particle Velocity 500 msec 900 msec

% close all
figure(2)
font_size = 20;

t = tiledlayout(1,4);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
t.Title.String = 'Particle Velocity at Free End';
t.Title.FontWeight = 'bold';
t.Title.FontSize = 28;


fig_size = [0.3 0.21 0.7 0.34];
nexttile
set(gcf,'color','w','units','normalized','position',fig_size,'MenuBar','none');
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
ylabel('Velocity [m/s]')
xlim([500 500.52])
ylim([0.015 0.042])
set(gca,'FontSize',font_size)

nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([600 600.52])
ylim([0.015 0.042])
set(gca,'FontSize',font_size)


nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3)
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3)
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([700 700.52])
ylim([0.015 0.042])
set(gca,'FontSize',font_size)

nexttile
plot(t_sim*10^3, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',3,'DisplayName','Model')
hold on
plot(t_data*10^3, v_data,':','Color','k','LineWidth',3,'DisplayName','Data')
hold off
grid on
% title('Particle Velocity at Free End')
xlabel('Time [msec]')
% ylabel('Velocity [m/s]','interpreter','latex')
xlim([800 800.52])
ylim([0.015 0.042])
set(gca,'FontSize',font_size)
legend

%% Long Term
close all
font_size = 24;


fig_size = [0.1 0.5 0.86 0.44];

figure(1)
set(gcf,'color','w','units','normalized','position',fig_size);
plot(t_sim, dudt_free*SCALE_RESULTS,'Color','r','LineWidth',0.3,'DisplayName','Model')
hold on
plot(t_data, v_data,':','Color','k','LineWidth',0.3,'DisplayName','Data')
hold off
grid on
title('Particle Velocity at Free End')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')
xlim([0 1])
set(gca,'FontSize',font_size)
legend

% figure(2)
% set(gcf,'color','w');
% set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44],'MenuBar','none');
% patchline(t_data, v_data,'edgecolor','k','linewidth',0.05,'edgealpha',0.05)
% hold on
% patchline(t_sim, dudt_free,'edgecolor','r','linewidth',0.05,'edgealpha',0.05)
% hold off
% title('Velocity vs. Time')
% xlabel('Time (seconds)','interpreter','latex')
% ylabel('Velocity (m/s)','interpreter','latex')
% xlim([0,1]) 
% grid on
% set(gca,'FontSize',font_size)

%% 3-D

t_sim_interp = linspace(0,1,length(t_sim));
v_right_interp = interp1(t_sim, dudt_free*SCALE_RESULTS, t_sim_interp);
[sim_FFT,  sim_freqs, sim_fs] =  MyFFT(v_right_interp,t_sim_interp); 
[data_FFT,data_freqs, data_fs] = MyFFT(v_data, t_data); 

close all
fig_size = [0.42 0.5 0.56 0.44];
figure(1)
set(gcf,'color','w','units','normalized','position',fig_size);%,'MenuBar','none');
plot(data_freqs*10^-3, data_FFT,'Color','k','linewidth',3)
hold on
plot(sim_freqs*10^-3, sim_FFT,'Color','r','linewidth',3)
hold off
title('FFT of Velocity at Free End')
xlabel('Frequency [kHz]','interpreter','latex')
ylabel('Amplitude','interpreter','latex')
xlim([0 300])
ylim([0.02*10^-6 2*10^-2])
grid on
legend('Model', 'Data')
set(gca,'FontSize',font_size,'YScale', 'log')

