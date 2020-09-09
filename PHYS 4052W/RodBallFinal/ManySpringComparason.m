clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')

t_data = t_post_impact; 
u_data = du_no_gravity; v_data = v_no_gravity;
c = v_phase;

%%
r = linspace(0,1,8)';
g = zeros(8,1);
b = linspace(1,0,8)';
rgb = [r,g,b];
font_size = 22;

close all
figure(1)
fig_size = [0.42 0.5 0.56 0.44];
set(gcf,'color','w','units','normalized','position',fig_size);

% 3 masses
N = 3;
kball = 13.5*10^9;
tspan = [0, 0.002];
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(1,:),'LineWidth',3,'DisplayName','N=3')
hold on
% 5 masses
N = 5;
kball = 16.9*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(2,:),'LineWidth',3,'DisplayName','N=5')
% 7 masses
N = 7;
kball = 16.9*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(3,:),'LineWidth',3,'DisplayName','N=7')
% 9 masses
N = 9;
kball = 16.5*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(4,:),'LineWidth',3,'DisplayName','N=9')
plot(t_data*10^3, u_data*10^3,'Color','g','LineWidth',7,'DisplayName','Data')
% 11 masses
N = 11;
kball = 16.2*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(5,:),'LineWidth',3,'DisplayName','N=11')
% 13 masses
N = 13;
kball = 16*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(6,:),'LineWidth',3,'DisplayName','N=13')

% 15 masses
N = 15;
kball = 16*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(7,:),'LineWidth',3,'DisplayName','N=15')

% 17 masses
N = 17;
kball = 15.7*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
plot(t, u,'Color',rgb(8,:),'LineWidth',3,'DisplayName','N=17')

hold off
title('Comparason of Mass-Spring Models')
xlabel('Time [seconds]')
ylabel('Displacement [mm]')
xlim([0 1])
ylim([-0.001 0.06])
grid on
set(gca,'FontSize',font_size)
lgd = legend;
lgd.Location = 'northwest';



%%

n_plots = 8;
r = linspace(0,1,n_plots)';
g = zeros(n_plots,1);
b = linspace(1,0,n_plots)';
rgb = [r,g,b];
font_size = 22;


f = @(v,N) v.*(N/9.6);

close all
figure(1)
fig_size = [0.42 0.5 0.56 0.44];
set(gcf,'color','w','units','normalized','position',fig_size);

% 3 masses
N = 3;
kball = 13.5*10^9;
tspan = [0, 0.002];
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(1,:),'LineWidth',3)
hold on
% 5 masses
N = 5;
kball = 16.9*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(2,:),'LineWidth',3)
% 7 masses
N = 7;
kball = 16.9*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(3,:),'LineWidth',3)
% 9 masses
N = 9;
kball = 16.5*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(4,:),'LineWidth',3)
plot(t_data*10^3, v_data,'Color','g','LineWidth',3)
% 11 masses
N = 11;
kball = 16.2*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(5,:),'LineWidth',3)
% 13 masses
N = 13;
kball = 16*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(6,:),'LineWidth',3)

% 15 masses
N = 15;
kball = 16*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(7,:),'LineWidth',3)

% 17 masses
N = 17;
kball = 15.7*10^9;
[t, u, v] = getSpringSol(N, kball, tspan);
v = f(v,N);
plot(t, v,'Color',rgb(8,:),'LineWidth',3)

hold off
title('Comparason of Mass-Spring Models')
xlabel('Time [seconds]')
ylabel('Displacement [mm]')
xlim([0 1])
ylim([-0.02 0.1])
grid on
set(gca,'FontSize',font_size)




