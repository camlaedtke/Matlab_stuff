clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')

%% Getting decay parameter from envelope

du = du_pendulum_subtracted;
% du = du_no_gravity;
t = t_post_impact;
fs = 1/(t(2)-t(1));
dudt = diff(du)./diff(t);
dudt(end+1) = 0;
% dudt = v_no_gravity;


dt = 1:497500;
sample_size = 1000;
[upper, lower] = envelope(du(dt),sample_size,'peak');
font_size = 28;

[fitobj, ~, ~] = fit(t(dt), upper, 'exp1');
fitobj
% Get envelope of signal
close all
figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44]);
% patchline(t(dt), du(dt),'edgecolor','r','linewidth',0.05,'edgealpha',0.05,'DisplayName','none')
plot(t(dt), du(dt),':','LineWidth',0.3,'Color', 'r')
hold on
p1 = plot(t(dt), upper, 'LineWidth', 3,'Color', 'b'); % upper envelope seems best
plot(t(dt), lower, 'LineWidth', 3,'Color', 'b')
p2 = plot(t(dt), fitobj(t(dt)), 'LineWidth', 4,'Color', 'g');
hold off
grid on
title('Displacement Data Pendulum Subtracted')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
ylim([-2.5*10^-6 2.5*10^-6])
set(gca,'FontSize', font_size)
legend([p1 p2],{'Envelope','Envelope Fit'})

%%

[pks, locs] = findpeaks(upper, t(dt), 'MinPeakDistance', tau/2);
[fitobj, ~, ~] = fit(locs, pks, 'exp1');
fitobj
figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
p1 = plot(locs, pks, 'LineWidth', 3,'Color', 'b'); % upper envelope seems best
hold on
% plot(locs, upper, 'LineWidth', 3,'Color', 'b') % upper envelope seems best
p2 = plot(locs, fitobj(locs), 'LineWidth', 3,'Color', 'g'); % upper envelope seems best
grid on
title('Displacement Response vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize', font_size)
legend([p1 p2],{'Envelope Peaks','Fit to Envelope Peaks'})


beta1 = 4.855;




%%

logdec = @(x) getdel(x);
D = logdec(pks);
D(end+1) = 0;

zeta = @(x) getzeta(x);
Z = zeta(pks);
Z(end+1) = 0;

delta_avg = mean(D);
zeta_avg = mean(Z);
fprintf('delta = %.5f, zeta = %.5f \n', delta_avg, zeta_avg)

close all
figure(2)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
plot(locs, D, 'LineWidth', 3,'Color', 'b') % upper envelope seems best
hold on
plot(locs, Z, 'LineWidth', 3,'Color', 'r') % upper envelope seems best
hold off
grid on
title('UPPER ENVELOPE: Displacement Response vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
set(gca,'FontSize', font_size)


%% % Logarithmic decrements: use peak amplitudes
close all
% USING ANALYTIC ENVELOPE
dt = 1:497500;
sample_size = 2000;
[upper, lower] = envelope(du(dt),sample_size,'analytic');
peak_dist = 0.0002;

% Get envelope of signal
figure(5)
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44],'MenuBar','none');
patchline(t(dt), du(dt),'edgecolor','r','linewidth',0.05,'edgealpha',0.05)
hold on
plot(t(dt), upper, 'LineWidth', 0.3,'Color', 'b')
plot(t(dt), lower, 'LineWidth', 0.3,'Color', 'b')
[pks, locs] = findpeaks(upper, t(dt));%,'MinPeakDistance',peak_dist);
findpeaks(upper, t(dt),'MinPeakDistance',peak_dist)
hold off
grid on
title('UPPER ENVELOPE: Displacement Response vs. Time')
xlabel('Time (seconds)')
% xlim([0.2 0.21])
ylabel('Displacement (m)')
set(gca,'FontSize', font_size)
x = pks;
[delta, zeta] = dampingratio(x);
figure(6)
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
plot(locs(1:end-1), delta,'LineWidth', 2,'Color', 'r','DisplayName', 'delta')
hold on
plot(locs(1:end-1), zeta,'LineWidth', 2,'Color', 'b','DisplayName', 'delta')
hold off
title('delta, zeta')
set(gca,'FontSize', font_size,'YScale','log')
legend
delta_avg = mean(delta);
zeta_avg = mean(zeta);
fprintf('Analytic, upper: delta = %.4e, zeta = %.4e \n', delta_avg, zeta_avg)





function [logdec, zeta] = dampingratio(x)
    logdec = zeros(length(x)-1,1);
    zeta = zeros(length(x)-1,1);
    for n = 1:length(x)-1
        x0 = x(n);
        x1 = x(n+1);
        logdec(n) = (1/n)*log(x0/x1);
        zeta(n) = 1/sqrt(1 + ((2*pi)/logdec(n))^2);
    end
end

function logdec = getdel(x)
    logdec = zeros(length(x)-1,1);
    for n = 1:length(x)-1
        x0 = x(n);
        x1 = x(n+1);
        logdec(n) = (1/n)*log(x0/x1);
    end
end

function zeta = getzeta(x)
    zeta = zeros(length(x)-1,1);
    for n = 1:length(x)-1
        x0 = x(n);
        x1 = x(n+1);
        logdec(n) = (1/n)*log(x0/x1);
        zeta(n) = 1/sqrt(1 + ((2*pi)/logdec(n))^2);
    end
end









