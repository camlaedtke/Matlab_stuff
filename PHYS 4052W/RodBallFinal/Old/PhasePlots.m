clear all
close all
load('tek0039PendulumRemoved.mat')
load('tek0039ImportantParams.mat')
load('tek0039Sim100MFreeEndData.mat')
%% Phase Space
close all
dt = 1:10:2400000;
font_size = 28;
view_shift = [160,10];

figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44]);
plot3(t_sim(dt), dudt_free(dt), dudx_free(dt),':','Color','r','LineWidth',0.8)
title('3-D Phase SPace')
xlabel('Time')
ylabel('Particle Velocity')
zlabel('Strain')
grid on
view(view_shift)
set(gca,'FontSize',font_size)


dt = 1:2:490000;
du = sgolayfilt(du_no_gravity,3,21);
dudt = diff(du)./diff(t_post_impact);
d2udt2 = diff(dudt)./diff(t_post_impact(1:end-1));

figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.01 0.56 0.44]);
plot3(t_post_impact(dt), dudt(dt), (A_rod/Y)*d2udt2(dt),':','Color','k','LineWidth',0.3)
title('3-D Phase SPace')
xlabel('Time')
ylabel('Velocity')
zlabel('Strain')
grid on
view(view_shift)
set(gca,'FontSize',font_size)


%% Cool stuff
close all
Fs_data = 1/(t_post_impact(2) - t_post_impact(1));

t_sim_interp = linspace(0,1,length(t_sim));
dudx_free_interp = interp1(t_sim, dudx_free, t_sim_interp);
Fs_model = 1/(t_sim_interp(2) - t_sim_interp(1));

figure(1)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44],'MenuBar','none');
pspectrum(d2udt2(dt),Fs_data,'persistence','FrequencyLimits',[0 400000],'TimeResolution',0.01)
set(gca,'FontSize',font_size)

figure(2)
set(gcf,'color','w');
set(gcf,'units','normalized','position',[0.42 0.05 0.56 0.44],'MenuBar','none');
pspectrum(dudx_free_interp,Fs_model,'persistence','FrequencyLimits',[0 400000],'TimeResolution',0.01)
set(gca,'FontSize',font_size)


