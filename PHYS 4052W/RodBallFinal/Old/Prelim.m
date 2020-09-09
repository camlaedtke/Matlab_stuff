%% Ball Rod Priliminary

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
vel_ave = diff(du_ave);%calculate the average displacement velocity

% GET V OF CENTER OF MASS
x = t_ave(dt);
y = du_ave(dt);
[fitobj, gof, outp] = fit(x, y,'poly1');

v_cm = fitobj.p1;
du_0 = du_ave(1);
du_ave = du_ave - du_0;

figure(1)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44]);
plot(t_ave(dt), du_ave(dt), 'LineWidth', 3 ,'Color', 'blue')
hold on
plot(t_ave(dt), fitobj(t_ave(dt)),'--','LineWidth', 4 ,'Color', 'black')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0 0.001]) 
ylim([0 0.00003])
grid on
legend('Data', 'Fit')
set(gca,'FontSize',font_size)

figure(2)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
plot(t_ave(dt), 10^6*vel_ave(dt), 'LineWidth', 3, 'Color', 'blue')
grid on
title('Velocity vs. Time')
xlabel('Time After Impact (seconds)')
ylabel('Velocity (m/sec)')
xlim([0 0.001]) 
ylim([-0.02 0.2])
set(gca,'FontSize',font_size)
[v_at_t_peaks, t_peaks] = findpeaks(10^6*vel_ave(dt), t_ave(dt), ...
    'MinPeakHeight', 0.02, 'MinPeakDistance',0.1*10^-3);
fprintf('\n t_peak = %.3e', t_peaks(1))
tau = mean(diff(t_peaks));
L = 0.5510;
fprintf('\n tau = %.3e\n', tau)

%%
close all

syms x(t) [128 1] 
syms omega_0

N_m = 128;
N_s = N_m-1;
M = x(t);
eqs = sym(ones(N_m,1));

eqs(1) = diff(M(1),t,2) == -omega_0^2*(M(1) - M(2));
for N = 2:(N_m-1)
    eqs(N) = diff(M(N),t,2) == -omega_0^2*(-M(N-1) + 2*M(N) - M(N+1));
end
eqs(N_m) = diff(M(N_m),t,2) == -omega_0^2*(M(N_m) - M(N_m-1));

[V] = odeToVectorField(eqs);
%%
w0 = N_m*sqrt(tau^2/(4*L^3))*10^8;
[V] = subs(V,omega_0, w0);
M = matlabFunction(V, 'vars', {'t','Y'});

v_0 = 6*v_cm;
Y_0 = zeros(1,2*N_m);
Y_0(2) = v_0;
tspan = [0, 0.01];
% tspan = t_ave(dt);

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
S = ode45(M, tspan, Y_0, options);
t_spr = S.x;
du_spr = S.y(2*N_m-1,:);
v_spr = S.y(2*N_m,:);

figure(3)
set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44]);
plot(t_ave(dt), du_ave(dt), 'LineWidth', 2, 'Color', 'blue')
hold on
plot(t_spr, du_spr, 'LineWidth', 2,'Color', 'red')
hold off
title('Displacement vs. Time')
xlabel('Time (seconds)')
ylabel('Displacement (m)')
xlim([0 0.002]) 
% ylim([-0.02 0.2])
grid on
set(gca,'FontSize',font_size)

figure(4)
set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44]);
plot(t_ave(dt), 10^6*vel_ave(dt), 'LineWidth', 2, 'Color', 'blue')
hold on
plot(t_spr, v_spr, 'LineWidth', 2,'Color', 'red')
hold off
title('Velocity vs. Time')
xlabel('Time (seconds)')
ylabel('Velocity (m/s)')
xlim([0 0.001]) 
% ylim([-0.02 0.2])
grid on
set(gca,'FontSize',font_size)

%%
save('tek0039PrelimParams.mat','t_ave','dt','du_ave','V_ave','vel_ave','v_cm','w0')
















