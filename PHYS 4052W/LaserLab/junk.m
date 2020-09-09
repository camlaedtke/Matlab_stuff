data_max = load('max_locations.mat', 'x_max');
data_min = load('min_locations.mat', 'x_min');

% Convert Units (Motor Steps) to Units (Meters)
params.steps_per_rev = 400; % steps/revolution
params.dist_per_rev = 0.00508; % meters/revolution
params.step_size = params.dist_per_rev/params.steps_per_rev; % meters/step

% Specify Measured Variables, Uncertainties

params.L_measured = 0.943 + 0.011; % meters
params.y_measured = 0.364 + 0.011; % meters

params.L_err = 0.004; % meters
params.y_err = 0.004; % meters
params.xM_err = 8.9*10^-6;% meters 
%  8.8e-6 for optimised r. chi squared, 
%  1e-5 for min, 7.2e-6 for max

params.F = sqrt((params.L_measured*(params.L_measured...
    -params.y_measured))/(2*params.y_measured));
params.xMax = data_max.x_max*params.step_size; % meters
params.xMin = data_min.x_min*params.step_size; % meters
params.mMax = linspace(0, length(params.xMax)-...
    1,length(params.xMax))';
params.mMin = linspace(0, length(params.xMin)...
    -1,length(params.xMin))';
params.GmMax = sqrt(4.*params.mMax+3/2);
params.GmMin = sqrt(4.*params.mMin+7/2);
params.xMax_err = ...
    params.xM_err.*ones(length(params.xMax),1); % meters
params.xMin_err = ...
    params.xM_err.*ones(length(params.xMin),1); % meters
FitPMTData(params, false)

x_vals = -0.01:0.0001:0.01;
for i = 1:length(x_vals)
    x = x_vals(i);
    params.L_measured = 0.943 + 0.011; % meters
    params.y_measured = 0.364 + 0.011; % meters
    params.L_measured = params.L_measured;
    params.y_measured = params.y_measured + x;
    FitPMTData(params, false)
   
end