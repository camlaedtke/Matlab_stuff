%% Import Data

clear all
data_max = load('max_locations.mat', 'x_max');
data_min = load('min_locations.mat', 'x_min');

%% Convert Units (Motor Steps) to Units (Meters)

steps_per_rev = 400; % steps/revolution
dist_per_rev = 0.00508; % meters/revolution
step_size = dist_per_rev/steps_per_rev; % meters/step
step_size*1000;

%% Specify Measured Variables, Uncertainties

L_measured = 0.943 + 0.011; % meters
y_measured = 0.364 + 0.011; % meters

F = sqrt((L_measured*(L_measured-y_measured))/(2*y_measured));
xMax = data_max.x_max*step_size; % meters
xMin = data_min.x_min*step_size; % meters
mMax = linspace(0, length(xMax)-1,length(xMax))';
mMin = linspace(0, length(xMin)-1,length(xMin))';
GmMax = sqrt(4.*mMax+3/2);
GmMin = sqrt(4.*mMin+7/2);

L_err = 0.002 + 0.002; % meters
y_err = 0.002 + 0.002; % meters
xM_err = 9*10^-6; % meters 8.8e-6 for both, 1e-5 for min,7.2e-6 for max
xMax_err = xM_err.*ones(length(xMax),1); % meters
xMin_err = xM_err.*ones(length(xMin),1); % meters
step_err = xM_err/step_size;
 
%% Perform LSQ Fit

max_fit = LSQfit(GmMax, xMax, xMax_err);
min_fit = LSQfit(GmMin, xMin, xMin_err);
Gm1 = GmMax.*max_fit(1,2)+max_fit(1,1);
Gm2 = GmMin.*min_fit(1,2)+min_fit(1,1);


%% Residuals

digits(4)
max_chivals = ChiVals(xMax,Gm1,xMax_err);
max_chi = max_chivals(:,1);
max_chi2 = max_chivals(:,2);
max_chi2_sum = sum(max_chi2);
max_dof = length(max_chi2) - 2;
max_chi2_reduced = vpa(max_chi2_sum/max_dof);

min_chivals = ChiVals(xMin,Gm2,xMin_err);
min_chi = min_chivals(:,1);
min_chi2 = min_chivals(:,2);
min_chi2_sum = sum(min_chi2);
min_dof = length(min_chi2) - 2;
min_chi2_reduced = vpa(min_chi2_sum/min_dof);

p_exceed_max = 1 - chi2cdf(max_chi2_sum,max_dof);
p_exceed_min = 1 - chi2cdf(min_chi2_sum,min_dof);

%% Plot Data

% clf
% fontsize = 20;
% f1 = get(groot, 'CurrentFigure');
% f1.Position = [100 100 800 800];
% e1 = errorbar(GmMax,xMax, xMax_err, 'o','DisplayName','Maxima');
% e1.Marker = 'o';
% e1.MarkerSize = 1;
% e1.Color = 'blue';
% e1.CapSize = 5;
% hold on
% plot(GmMax,Gm1,'b','DisplayName','Maxima Fit')
% e2 = errorbar(GmMin,xMin, xMin_err, 's','DisplayName','Minima');
% e2.Marker = 'o';
% e2.MarkerSize = 1;
% e2.Color = 'red';
% e2.CapSize = 5;
% plot(GmMin,Gm2,'r','DisplayName','Minima Fit')
% title('Plot of Peak Distances $X_{Maxima}$, $X_{Minima}$ vs. $G(m, k)$',...
%     'FontSize',fontsize,'Interpreter','latex')
% xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
% ylabel('Distance (meters)','FontSize',fontsize,'Interpreter','latex')
% legend('Location','northwest','FontSize',fontsize,'Interpreter','latex',...
%     'NumColumns',2)
% hold off
% 
% f2 = figure;
% hold on
% scatter(GmMax, max_chi,'filled');
% ylim([-3 3])
% f2.Position = [900 100 600 400];
% title('Weighted Residuas for Maxima Fit',...
%     'FontSize',fontsize,'Interpreter','latex')
% xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
% ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
% grid on
% hold off
% 
% f3 = figure;
% hold on
% scatter(GmMin,min_chi,'filled','r');
% ylim([-3 3])
% f3.Position = [900 500 600 400];
% title('Weighted Residuals for Minima Fit',...
%     'FontSize',fontsize,'Interpreter','latex')
% xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
% ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
% grid on
% hold off
% 
% hold on
% scatter(GmMin,min_chi,'filled','r');
% ylim([-3 3])
% f3.Position = [900 500 600 400];
% title('Weighted Residuals for Minima Fit',...
%     'FontSize',fontsize,'Interpreter','latex')
% xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
% ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
% grid on
% hold off
% 
% f4 = figure;
% hold on
% scatter(GmMin,min_chi2,'filled','r');
% scatter(GmMax,max_chi2,'filled','b');
% f4.Position = [1500 100 800 600];
% title('Weighted $\chi^2$ values for Both Fits',...
%     'FontSize',fontsize,'Interpreter','latex')
% xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
% ylabel('$\chi^2$ value','FontSize',fontsize,'Interpreter','latex')
% grid on
% hold off

%% Derive Error Propagation Equation for Lambda

syms F G L y m k lambda x_m x_0 s
assume(y>0)
assume(L>0)
assume(m>0)
assume(x_m>0)
assume(F>0)
assume(G>0)
assume(s>0)
assume(lambda>0)
x_m = sqrt(lambda)*F*G;
x_m = x_m+x_0; % MAIN EQUATION
equ = x_m == subs(x_m, sqrt(lambda)*F, s);
x_m = subs(x_m, sqrt(lambda)*F, s);
lambda = solve(equ, lambda);
lambda = subs(lambda,F,sqrt((L*(L-y))/(2*y)));
syms Dlambdads sigma_s DlambdaDL sigma_L DlambdaDy sigma_y
DlambdaDs = simplify(diff(lambda, s));
DlambdaDL = simplify(diff(lambda, L));
DlambdaDy = simplify(diff(lambda, y));
sigma_lambda = sqrt(sigma_s^2*(DlambdaDs)^2+...
                    sigma_L^2*(DlambdaDL)^2+...
                    sigma_y^2*(DlambdaDy)^2);
                
%% Calculate Lambda (Maxima, Minima)

digits(3)  
s_max = max_fit(1,2);
s_max_err = max_fit(2,2);
old = [s, L, y];
new = [s_max, L_measured, y_measured];
fprintf('Maxima Fit Slope Uncertainty: %.2e \n',s_max)
lambda_max_data = vpa(subs(lambda, old, new));
s_min = min_fit(1,2);
s_min_err = min_fit(2,2);
new = [s_min, L_measured, y_measured];
lambda_min_data = vpa(subs(lambda, old, new));

% Calculate Error in Lambda (maxima, minima)
old = [s, L, y, sigma_s, sigma_L, sigma_y];
new = [s_max, L_measured, y_measured, s_max_err, L_err, y_err];
fprintf('Minima Fit Slope Uncertainty: %.2e \n',s_min)
lambda_max_err = vpa(subs(sigma_lambda, old, new));
new = [s_min, L_measured, y_measured, s_min_err, L_err, y_err];
lambda_min_err = vpa(subs(sigma_lambda, old, new));

% Weighted Average
lambda_obs = [lambda_max_data, lambda_min_data];
lambda_obs_errs = [lambda_max_err, lambda_min_err];
sum1 = sum(lambda_obs./lambda_obs_errs.^2);
sum2 = sum(1./lambda_obs_errs.^2);
lambda_weighted_avg = sum1/sum2;
lambda_weighted_avg_err = sqrt(1/sum2);

% Display Results
str1 = '+/-';
str2 = 'nm';
str3 = 'mm';
str4 = 'motor steps';
chi_str =  '        Chi Squared =';
chi_str2 = 'Reduced Chi Squared =';
p_val_str = '            P value =';

formatSpec = "Lambda From Maxima  = %g %s %.2g %s";

X1 = sprintf(formatSpec,lambda_max_data*10^9 ,str1,...
    lambda_max_err*10^9,str2);

formatSpec = "%s %.2f";
X1_chi = sprintf(formatSpec,chi_str ,max_chi2_sum);
X1_chi2 = sprintf(formatSpec,chi_str2,max_chi2_reduced);
X1_p_val = sprintf(formatSpec,p_val_str,p_exceed_max);

formatSpec = "Lambda From Minima  = %g %s %.2g %s";
X2 = sprintf(formatSpec,lambda_min_data*10^9 ,str1,...
    lambda_min_err*10^9,str2);

formatSpec = "%s %.2f";
X2_chi = sprintf(formatSpec,chi_str ,min_chi2_sum);
X2_chi2 = sprintf(formatSpec,chi_str2,min_chi2_reduced);
X2_p_val = sprintf(formatSpec,p_val_str,p_exceed_min);

lambda_nm = vpa(lambda_weighted_avg*10^9);
lambda_err_nm = vpa(lambda_weighted_avg_err*10^9);
formatSpec = "Lambda Weighted Avg = %.0f %s %.0f %s";
X3 = sprintf(formatSpec, lambda_nm,str1,lambda_err_nm,str2);

formatSpec = "Motor Step Distance: %.3f %s ";
X4 = sprintf(formatSpec,step_size*10^3,str3);

formatSpec = "Motor Step Error: %.3f %s or %.2g %s";
X5 = sprintf(formatSpec, xM_err*10^3, str3, step_err, str4);

formatSpec = "Measured y, L: %.4g , %.4g %s ";
X6 = sprintf(formatSpec,y_measured,L_measured,str3);

formatSpec = "Measured y, L Error: %.0f %s ";
X7 = sprintf(formatSpec,y_err*10^3,str3);

results = [X1;X1_chi;X1_chi2;X1_p_val;X2;X2_chi;X2_chi2;X2_p_val;X3];
params_and_errors = [X4;X5;X6;X7];
disp(results)
disp(params_and_errors)

