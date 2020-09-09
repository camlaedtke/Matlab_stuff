clear all
close all
% CALIBRATION DATA
attenuation = 1/998.2;
calibration_data = dlmread('Noise Lab.csv',',',1,0);
freq = calibration_data(:,1);
V_in = calibration_data(:,2);
V_in_err = calibration_data(:,3);
V_out = calibration_data(:,4);
V_out_err = calibration_data(:,5);
V_in = V_in.*attenuation;
V_in_err = V_in_err.*attenuation;
gain = V_out./V_in;

T = 295; % degrees Kelvin
T_err = 3; % degrees Kelvin

% ERROR IN TRAPIZOIDAL INTEGRATION
G2B_err = zeros(length(V_in)-1,1);
for i = 1:length(V_in)-1
    G2B_err(i) = sqrt((V_out_err(i)^2/V_in(i)^2 + ...
        (V_out(i)^2*V_in_err(i)^2)/(V_in(i)^4))*(freq(i)/2 - ...
        freq(i+1)/2)^2 + (V_out_err(i+1)^2/V_in(i+1)^2 + ...
        (V_out(i+1)^2*V_in_err(i+1)^2)/(V_in(i+1)^4))*(freq(i)/2 -...
        freq(i+1)/2)^2);
end
G2B_err = sum(G2B_err);

R_100k = 100000; % ohms
shot_100k = dlmread('ShotNoise.csv',',',[1 0 7 3]).*(1/1000);
V_ave_100k = shot_100k(:,1);
V_ave_err_100k = shot_100k(:,2);
V_rms_100k = shot_100k(:,3);
V_rms_err_100k = shot_100k(:,4);

R_200k = 200000; % ohms
shot_200k = dlmread('ShotNoise.csv',',',[10,0,16,3]).*(1/1000);
V_ave_200k = shot_200k(:,1);
V_ave_err_200k = shot_200k(:,2);
V_rms_200k = shot_200k(:,3);
V_rms_err_200k = shot_200k(:,4);

%% 100 K RESISTOR 
close all
G2B = trapz(freq,(gain.^2));

%x = V_rms_100k.^2;
%y = 2*R_100k.*V_ave_100k.*G2B;
y = V_rms_100k.^2;
x = 2*R_100k.*V_ave_100k.*G2B;

dydVave = 2*R_100k.*G2B;
dydG2B = 2*R_100k.*V_ave_100k;

%sig_y = sqrt(V_ave_err_100k.^2.*dydVave.^2+G2B_err.^2.*dydG2B.^2);
sig_y = V_rms_err_100k;


[fitobj, gof, outp] = fit(x,y,'poly1','Weights',(1./sig_y).^2);
fprintf('\nChi Squared R = 100k: %.5f', gof.sse);

figure(1)
set(gcf,'units','normalized','position',[0.5 0.5 0.4 0.5]);
y_plot1 = x.*fitobj.p1 + fitobj.p2;
%e = 1/fitobj.p1;
e = fitobj.p1;
plot(x,y_plot1)
hold on
title('Shot Noise Linear Fit with $ R = 100 k \Omega $','Interpreter','latex')
xlabel('$2RV_{ave} \int_{0}^{\infty} g^2(f) df$','Interpreter','latex')
ylabel('$V_{RMS}^2$','Interpreter','latex')
set(gca,'FontSize',22)
grid on
errorbar(x,y,sig_y,'o')
hold off

figure(2)
set(gcf,'units','normalized','position',[0.7 0.05 0.4 0.5]); 
plot(x(1:length(outp.residuals),:), outp.residuals, 'x','MarkerEdgeColor','b','MarkerSize',10)
title('Weighted Residuals','Interpreter','latex')
xlabel('$2RV_{ave} \int_{0}^{\infty} g^2(f) df$','Interpreter','latex')
ylabel('$\chi^2$','Interpreter','latex')
set(gca,'FontSize',20)
grid on


fprintf('\n100k Ohms e = %.3e pm %.3e \n',e)

mycoeffnames = coeffnames(fitobj);   
mycoeffvalues = coeffvalues(fitobj);
mycoeff_onesig = diff(confint(fitobj))/4; 
error_matrix = inv(outp.Jacobian'*outp.Jacobian); 
dcov =  diag(error_matrix);           
unc = sqrt(dcov)';             
disp( 'Fit Parameters and Uncertainties (+/- One Standard Deviation)')
for irow = 1:outp.numparam 
    fprintf('%s: %.4e +/- %.4e\n', char( mycoeffnames(irow)),...
        mycoeffvalues(irow), unc(irow)) 	%Jacobian Method
end

% DETERMINE THE AMPLIFIER NOISE AT THE INPUT

%% 200 K RESISTOR
close all
G2B = trapz(freq,(gain.^2));


y = V_rms_200k.^2;
x = 2*R_200k.*V_ave_200k.*G2B;

% dydVave = 2*R_200k.*G2B;
% dydG2B = 2*R_200k.*V_ave_200k;

sig_y = V_rms_err_200k;
              
[fitobj, gof, outp] = fit(x,y,'poly1','Weights',(1./sig_y).^2);
fprintf('\nChi Squared R = 200k: %.5f \n', gof.sse);
e = fitobj.p1;
fprintf('200k Ohms e = %.3e pm %.3e \n',e,mean(e_err_200k))


figure(2)
set(gcf,'units','normalized','position',[0.5 0.1 0.3 0.4]);
y_plot2 = x.*fitobj.p1 + fitobj.p2;
plot(x,y_plot2)
hold on
title('Shot Noise Linear Fit with $ R = 200 k \Omega $','Interpreter','latex')
xlabel('$2RV_{ave} \int_{0}^{\infty} g^2(f) df$','Interpreter','latex')
ylabel('$V_{RMS}^2$','Interpreter','latex')
set(gca,'FontSize',22)
grid on
errorbar(x,y,sig_y,'o')
hold off


figure(3)
set(gcf,'units','normalized','position',[0.7 0.05 0.4 0.5]); 
plot(x(1:length(outp.residuals),:), outp.residuals, 'x','MarkerEdgeColor','b','MarkerSize',10)
title('Weighted Residuals','Interpreter','latex')
xlabel('$2RV_{ave} \int_{0}^{\infty} g^2(f) df$','Interpreter','latex')
ylabel('$\chi^2$','Interpreter','latex')
set(gca,'FontSize',20)
grid on

mycoeffnames = coeffnames(fitobj);   
mycoeffvalues = coeffvalues(fitobj);
mycoeff_onesig = diff(confint(fitobj))/4; 
error_matrix = inv(outp.Jacobian'*outp.Jacobian); 
dcov =  diag(error_matrix);           
unc = sqrt(dcov)';             
disp( 'Fit Parameters and Uncertainties (+/- One Standard Deviation)')
for irow = 1:outp.numparam 
    fprintf('%s: %.4e +/- %.4e\n', char( mycoeffnames(irow)),...
        mycoeffvalues(irow), unc(irow)) 	%Jacobian Method
end
