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
length(gain)

figure(1)
set(gcf,'units','normalized','position',[0.7 0.5 0.4 0.5]);
scatter(freq,gain)
xlabel("Frequency (Hz)",'Interpreter','latex')
ylabel("Gain ($V_{out}/V_{in}$)",'Interpreter','latex')
title("Calibration Readings from 100 Hz to 40,000 Hz",'Interpreter','latex')
set(gca,'FontSize',22)
grid on

% LOW RESISTANCE MEASUREMENTS
JohnsonData = dlmread('JohnsonData.csv',',',1,0);
R_L = JohnsonData(:,1); % L stands for 'low' resistance
V_R_L = JohnsonData(:,2); 
V_R_L_err = JohnsonData(:,3);
V_L_zero = JohnsonData(1,2);
V_L_zero_err = JohnsonData(1,3);

% HIGH RESISTANCE MEASUREMENTS
Johnson2 = dlmread('JohnsonData2.csv',',',1,0);
R_H = Johnson2(:,1);
V_R_H = Johnson2(:,2);
V_R_H_err = Johnson2(:,3);
V_H_zero = Johnson2(:,4); 
V_H_zero_err = Johnson2(:,5);

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

%% LOW RESISTANCE - NO CORRECTION FACTOR
G2B = trapz(freq,(gain.^2));

x = R_L;
y = V_R_L.^2./(4*G2B*T);
dydVR = V_R_L./(2*G2B*T);
dydT = -(V_R_L.^2./(4*G2B*T^2));
dydG2B = -(V_R_L.^2./(4*G2B.^2*T));
sig_y = sqrt(V_R_L_err.^2.*dydVR.^2+T_err^2.*dydT.^2+G2B_err^2.*dydG2B.^2);

[fitobj, gof, outp] = fit(x,y,'poly1','Weights',(1./sig_y).^2);
fprintf('\nChi Squared Low R fit: %.5f', gof.sse);

figure(1) 
set(gcf,'units','normalized','position',[0.6 0.5 0.3 0.4]);
y_plot1 = x.*fitobj.p1 + fitobj.p2;
k_L = fitobj.p1;
plot(x,y_plot1);
hold on
title('Johnson Noise Linear Fit Low R','Interpreter','latex');
xlabel('Resistance (Ohms)','Interpreter','latex')
ylabel('$\frac{V_{R}^2}{4 T G^2 B}$','Interpreter','latex')
set(gca,'FontSize',28)
grid on
errorbar(x,y,sig_y,'o');
hold off

figure(2)
set(gcf,'units','normalized','position',[0.6 0.05 0.3 0.4]); 
plot(x(1:length(outp.residuals),:), outp.residuals, 'x','MarkerEdgeColor','b','MarkerSize',10)
title('Weighted Residuals vs. Resistance','Interpreter','latex')
xlabel('Resistance (Ohms)','Interpreter','latex')
ylabel('$\chi^2$','Interpreter','latex')
set(gca,'FontSize',28)
grid on

fprintf('\nLow R Measurement of k with NO correction factor: %.6g \n',k_L)

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

%% HIGH RESISTANCE VR MEASUREMENTS
close all 	
% CORRECTION FACTOR FOR STRAY CAPACITANCE
C0 = 30.39*10^-12; % Gives vest chi value
G2B = zeros(length(R_H),1);
for i = 1:length(R_H)
    G2B(i) = trapz(freq,(gain.^2)./(1+((2*pi.*freq.*R_H(i)*C0).^2)));
    %fprintf('G2B (R = %.0f Ohms) has %.7f %% Error \n',R_H(i),(G2B_err/G2B(i))*100)
end

x = R_H;
y = V_R_H.^2./(4*G2B*T);
dydVR = V_R_H./(2*G2B*T);
dydT = -(V_R_H.^2./(4*G2B*T^2));
dydG2B = -(V_R_H.^2./(4*G2B.^2*T));
sig_y = sqrt(V_R_H_err.^2.*dydVR.^2+T_err^2.*dydT.^2+G2B_err^2.*dydG2B.^2);

[fitobj, gof, outp] = fit(x, y,'poly1', 'Weights',(1./sig_y).^2);
fprintf('\nChi Squared High R fit: %.5f', gof.sse);

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


figure(3) 
set(gcf,'units','normalized','position',[0.7 0.5 0.3 0.4]); 
y_plot2 = x.*fitobj.p1 + fitobj.p2;
k_H = fitobj.p1;
plot(x,y_plot2);
hold on
title('Johnson Noise Fit With Correction Factor (High R)','Interpreter','latex');
xlabel('Resistance (Ohms)','Interpreter','latex')
ylabel('$\frac{V_{R}^2}{4 T G^2 B}$','Interpreter','latex')
set(gca,'FontSize',28)
grid on
errorbar(x,y,sig_y,'o');
hold off

figure(4)
set(gcf,'units','normalized','position',[0.7 0.05 0.3 0.4]); 
plot(x(1:length(outp.residuals),:), outp.residuals, 'x','MarkerEdgeColor','b','MarkerSize',10)
title('Weighted Residuals vs. Resistance','Interpreter','latex')
xlabel('Resistance (Ohms)','Interpreter','latex')
ylabel('$\chi^2$','Interpreter','latex')
set(gca,'FontSize',28)
grid on
fprintf('\nHigh R Measurement of k with C0 = %.3e correction factor: %.6g \n',C0,k_H)

%% GETTING THE UNCERTAINTY IN C0 BY MINIMIZING LEAST SQUARES
close all 	

% Range of C0 values
%a = 28.69:0.01:32.19;
a = 10:0.1:60;
C = (10^-12.*a)';



sse = zeros(length(C),1);



for i = 1:length(C)
    C0 = C(i);
    
    sse(i) = GetSSE(C0,R_H,freq,gain,V_R_H,V_R_H_err,T,T_err,G2B_err);
        
end

C_sse = [C,sse];
figure(5)
set(gcf,'units','normalized','position',[0.5 0.3 0.5 0.6]); 
scatter(C_sse(:,1),C_sse(:,2),5)
title('Correction Factor vs. Sum of Squares','Interpreter','latex')
xlabel('Correction Factor $C_{0}$','Interpreter','latex')
ylabel('$\Sigma \chi^2$','Interpreter','latex')
set(gca,'FontSize',28)
grid on

C_sse(:,2) = -1.*C_sse(:,2);
x = C_sse(:,1);
y = C_sse(:,2);

[px,py] = findpeaks(y,x);
findpeaks(y,x)


fprintf('\n\nBest Value of C0 is: C0 = %.4e\n',py)
fprintf('Minimum Sum of Chi Squared + 1 = %.4f \n',-px+1) % C(px+1) = 3.216 * 10^-11
C_1sigma = 3.216e-11 - py;
fprintf('Uncertainty in C0 is: C_err = %.3g \n',C_1sigma)


function Chival = GetSSE(C0,R_H,freq,gain,V_R_H,V_R_H_err,T,T_err,G2B_err)
    G2B = zeros(length(R_H),1);
    
    for i = 1:length(R_H)
        G2B(i) = trapz(freq,(gain.^2)./(1+((2*pi.*freq.*R_H(i)*C0).^2)));
    end

    x = R_H;
    y = V_R_H.^2./(4*G2B*T);
    dydVR = V_R_H./(2*G2B*T);
    dydT = -(V_R_H.^2./(4*G2B*T^2));
    dydG2B = -(V_R_H.^2./(4*G2B.^2*T));
    sig_y = sqrt(V_R_H_err.^2.*dydVR.^2+T_err^2.*dydT.^2+G2B_err^2.*dydG2B.^2);

    [~, gof, ~] = fit(x, y,'poly1', 'Weights',(1./sig_y).^2);
    fprintf('\nChi Squared Fit --- C = %.4e --- SSE = %.5f',C0,gof.sse);
    Chival = gof.sse;
end


