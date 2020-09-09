clear all;
close all
load('tek0039Averaged.mat')
font_size = 20;

t_impact = 0.00487;
ix_impact = find(abs(t_ave(dt)-t_impact) < 0.000001);
fprintf('\nIndex at t_impact: [ %g ]', ix_impact);
ix_impact = ix_impact (1);
t_ave = t_ave - t_impact;
dt = ix_impact:length(dt);
fprintf('\nAfter t_impact: %.0f data points \n',length(t_ave(dt)))


%%
vel_ave = diff(du_ave);%calculate the average displacement velocity
close all

figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44]);
plot(t_ave(dt), du_ave(dt), 'LineWidth', 3 ,'Color', 'blue')
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0 0.0011]) 
grid on
set(gca,'FontSize',font_size)

order = 3;
framelen = 11;
vel_ave_filt = sgolayfilt(vel_ave(dt),order, framelen);
fprintf('\nFiltered vel_ave: %.0f data points \n',length(vel_ave_filt))

figure(2)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
plot(t_ave(dt), vel_ave(dt), ':' ,'LineWidth', 3)
hold on
plot(t_ave(dt), vel_ave_filt,'.-', 'LineWidth', 3)
hold off
grid on
title('Velocity vs. Time')
xlabel('Time After Impact (seconds)')
ylabel('Velocity (m/sec)')
xlim([0.00023 0.00035]) 
legend('original','sgolay')
set(gca,'FontSize',font_size)