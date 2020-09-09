function FitPMTData(params, plot)
    % Perform LSQ Fit
    params.max_fit = LSQfit(params.GmMax, params.xMax, params.xMax_err);
    params.min_fit = LSQfit(params.GmMin, params.xMin, params.xMin_err);
    params.Gm1 = params.GmMax.*params.max_fit(1,2)+params.max_fit(1,1);
    params.Gm2 = params.GmMin.*params.min_fit(1,2)+params.min_fit(1,1);

    % Residuals
    digits(4)
    params.max_chivals = ChiVals(params.xMax,params.Gm1,params.xMax_err);
    params.max_chi = params.max_chivals(:,1);
    params.max_chi2 = params.max_chivals(:,2);
    params.max_chi2_sum = vpa(sum(params.max_chi2));
    params.max_dof = length(params.max_chi2) - 2;
    params.max_chi2_reduced = vpa(params.max_chi2_sum/params.max_dof);
    
    params.min_chivals = ChiVals(params.xMin,params.Gm2,params.xMin_err);
    params.min_chi = params.min_chivals(:,1);
    params.min_chi2 = params.min_chivals(:,2);
    params.min_chi2_sum = vpa(sum(params.min_chi2));
    params.min_dof = length(params.min_chi2) - 2;
    params.min_chi2_reduced = vpa(params.min_chi2_sum/params.min_dof);
    
    if plot == true
        PlotPMTData(params)
    else 
       
    end
    
    % Derive Error Propagation Equation for Lambda
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
    
    % Calculate Lambda (Maxima, Minima)

    digits(3)
    s_max = params.max_fit(1,2);
    s_max_err = params.max_fit(2,2);
    old = [s, L, y];
    new = [s_max, params.L_measured, params.y_measured];
    lambda_max_data = vpa(subs(lambda, old, new));
    
    s_min = params.min_fit(1,2);
    s_min_err = params.min_fit(2,2);
    new = [s_min, params.L_measured, params.y_measured];
    lambda_min_data = vpa(subs(lambda, old, new));
    
    % Calculate Error in Lambda (maxima, minima)
    old = [s, L, y, sigma_s, sigma_L, sigma_y];
    new = [s_max, params.L_measured, params.y_measured, s_max_err,...
        params.L_err, params.y_err];
    lambda_max_err = vpa(subs(sigma_lambda, old, new));
    
    new = [s_min, params.L_measured, params.y_measured, s_min_err,...
        params.L_err, params.y_err];
    lambda_min_err = vpa(subs(sigma_lambda, old, new));
    
    % Weighted Average
    lambda_obs = [lambda_max_data, lambda_min_data];
    lambda_obs_errs = [lambda_max_err, lambda_min_err];
    sum1 = sum(lambda_obs./lambda_obs_errs.^2);
    sum2 = sum(1./lambda_obs_errs.^2);
    params.lambda_weighted_avg = sum1/sum2;
    params.lambda_weighted_avg_err = sqrt(1/sum2);
    
    % Display Results
    
    
    fprintf('y = %.4g L = %.4g \n',params.y_measured, params.L_measured)
%     fprintf("Maxima = %.0f +/- %.0f nm \n",lambda_max_data*10^9 ,...
%         lambda_max_err*10^9);
% 
%     fprintf("Minima = %g +/- %.2g nm \n",lambda_min_data*10^9,...
%         lambda_min_err*10^9);
    
    lambda_nm = vpa(params.lambda_weighted_avg*10^9);
    lambda_err_nm = vpa(params.lambda_weighted_avg_err*10^9);
    fprintf("Lambda Avg = %.0f +/- %.0f nm \n\n", lambda_nm,lambda_err_nm);
     
  

end

