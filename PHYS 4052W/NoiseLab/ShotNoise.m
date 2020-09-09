% Set the switch on shielded box to "shot" position
% photodiode connected across resistor and BNC output
% Adjusted LED intensity and observe the range of Vave 
% as a DC voltage across R

% Range of Vave: 0.00001V to 0.185 V

R_100k = 100000; % ohms
shot_100k = dlmread('ShotNoise.csv',',',[1 0 7 3]).*(1/1000);
V_ave_100k = shot_100k(:,1);
V_ave_err_100k = shot_100k(:,2);
V_rms_100k = shot_100k(:,3);
V_rms_Err_100k = shot_100k(:,4);

x = 2*R_100k.*V_ave_100k.*intg2C;
y = V_rms_100k.^2;
scatter(x,y)

xlabel('Resistance (Ohms)')
ylabel('$\frac{V_{R}^2}{4 T G^2 B}$','Interpreter','latex')
fit_vals = LSQfit(x,y,V_rms_Err_100k);
slope1 = fit_vals(1,2);
intercept1 = fit_vals(1,1);

%%
R_200k = 200000; % ohms
shot_200k = dlmread('ShotNoise.csv',',',[10,0,16,3]).*(1/1000);
V_ave_200k = shot_200k(:,1);
V_ave_err_200k = shot_200k(:,2);
V_rms_200k = shot_200k(:,3);
V_rms_Err_200k = shot_200k(:,4);

x = 2*R_200k.*V_ave_200k.*intg2C;
y = V_rms_200k.^2;
scatter(x,y)

fit_vals = LSQfit(x,y,V_rms_Err_200k);
slope2 = fit_vals(1,2)
intercept2 = fit_vals(1,1)