%% Mossbauer Lab

% Conversion factor
p1 = 450.06928;
p1_err = 0.0853;
p6 = 578.1798;
p6_err = 0.11694;
delta_p = p6-p1;
delta_p_err = sqrt(p1_err^2 + p6_err^2);

delta_v = 10.657;
delta_v_err = 0.017;
C = delta_v/delta_p;
C_err = sqrt(delta_v_err^2*(1/delta_p)^2 + delta_p_err^2*(delta_v/delta_p^2)^2);
offset = - (p6 - delta_p/2);
offset_err = sqrt(delta_p_err^2*(1/2)^2 + p6_err^2);
C_err = C_err + offset_err;

% GET DELTA 0
% 1) Get P4, P2 in terms of V4, V2
p4 = 524.07792;
p4_err = 0.15999;
p2 = 477.23072;
p2_err = 0.09698;

% 2) Put P4, P2 into sigma_C equation to get V4_err, V2_err
v4 = (p4+offset)*C;
v4_err = p4_err*C_err;
v2 = (p2+offset)*C;
v2_err = p2_err*C_err;

% 3) Put V4, V4_err into sigma_E equation to get E4_err
E_gamma = 14410; %eV
E_gamma_err = 20;
c = 299792458;

E4 = (E_gamma/(1000*c))*v4;
E4_err = (1/(1000*c))*sqrt(E_gamma_err^2*v4^2 + v4_err^2*E_gamma^2);
E2 = (E_gamma/(1000*c))*v2;
E2_err = (1/(1000*c))*sqrt(E_gamma_err^2*v2^2 + v2_err^2*E_gamma^2);

delta_0 = E4-E2;
delta_0_err = sqrt(E4_err^2+E2_err^2);

% GET DELTA 1
% 1)
p6 = 578.1798;
p6_err = 0.11694;
p5 = 551.15342;
p5_err = 0.12363;

% 2)
v6 = (p6+offset)*C;
v6_err = p6_err*C_err;
v5 = (p5+offset)*C;
v5_err = p5_err*C_err;

% 3)
E6 = (E_gamma/(1000*c))*v6;
E6_err = (1/(1000*c))*sqrt(E_gamma_err^2*v6^2 + v6_err^2*E_gamma^2);
E5 = (E_gamma/(1000*c))*v5;
E5_err = (1/(1000*c))*sqrt(E_gamma_err^2*v5^2 + v5_err^2*E_gamma^2);

delta_1 = E6-E5;
delta_1_err = sqrt(E6_err^2+E5_err^2);
fprintf(['\n------ CALCULATIONS -------\n\n'....
         'Delta_1 = %.3e +/- %.1e eV \n'...
         'Delta_0 = %.3e +/- %.1e eV'],delta_1,delta_1_err,...
    delta_0,delta_0_err)


% RATIO OF MAGNETIC MOMENTS
moment_ratio = delta_0/delta_1;
moment_ratio_err = sqrt(delta_0_err^2*(1/delta_1)^2+...
    delta_1_err^2*(delta_0/delta_1^2)^2);
fprintf('\nmu_0/mu_1 = %.3f +/- %.3f', moment_ratio, moment_ratio_err)

% SOLVE FOR MAGNETIC FIELD AT NUCLEUS
mu_N = 3.152451*10^-8; % eV
mu_0 = 0.0903*mu_N;
mu_0_err = 0.0007*mu_N;
mu_1 = 0.1548*mu_N;
mu_1_err = 0.0013*mu_N;

H0 = delta_0/(2*mu_0);
H0_err = sqrt(delta_0_err^2*(1/(2*mu_0))^2+...
                mu_0_err^2*(delta_0/(2*mu_0^2))^2);
H1 = (3*delta_1)/(2*mu_1);
H1_err = sqrt(delta_1_err^2*(3/(2*mu_1))^2+...
            mu_1_err^2*((3*delta_1)/(2*mu_1^2))^2);
fprintf(['\nH_0 = %.2f +/- %.2f Tesla \n'...
          'H_1 = %.2f +/- %.2f Tesla \n'],H0,H0_err,H1,H1_err)

% Preliminary vals
d_1_accepted = 1.08e-7;
d_0_accepted = 1.88e-7;
moment_accepted = 1.75;
H_accepted = 33;

s_delta_1 = get_sigma(delta_1, d_1_accepted, delta_1_err);
s_delta_0 = get_sigma(delta_0, d_0_accepted, delta_0_err);
s_moment_ratio = get_sigma(moment_ratio, moment_accepted, moment_ratio_err);
s_H_1 = get_sigma(H1, H_accepted, H1_err);
s_H_0 = get_sigma(H0, H_accepted, H0_err);

fprintf(['\n    sigma Delta_1 = %.2f \n'...
           '    sigma Delta_0 = %.2f \n'...
           '    sigma mu_0/mu_1 = %.2f \n'...
           '    sigma_H_0 = %.2f \n'...
           '    sigma_H_1 = %.2f \n\n'],...
            s_delta_1,...
            s_delta_0,...
            s_moment_ratio,...
            s_H_1,...
            s_H_0)
         
% IMPROVE ON VALUE OF DELTA_0 AND DELTA_1
% Get remaining energies
% 1)
p1 = 450.06928;
p1_err = 0.0853;
p3 = 503.83584;
p3_err = 0.17905;

% 2)
v1 = (p1+offset)*C;
v1_err = p1_err*C_err;
v3 = (p3+offset)*C;
v3_err = p3_err*C_err;

% 3)
E1 = (E_gamma/(1000*c))*v1;
E1_err = (1/(1000*c))*sqrt(E_gamma_err^2*v1^2 + v1_err^2*E_gamma^2);
E3 = (E_gamma/(1000*c))*v3;
E3_err = (1/(1000*c))*sqrt(E_gamma_err^2*v3^2 + v3_err^2*E_gamma^2);

E = [E1, E2, E3, E4, E5, E6]';

A = [3/2,  -1/2,  1,   1;
     1/2,  -1/2,  1,  -1;
    -1/2,  -1/2,  1,  -1;
     1/2,   1/2,  1,  -1;
    -1/2,   1/2,  1,  -1;
    -3/2,   1/2,  1,   1];
    
vars = A\E;
delta_1 = -vars(1);
delta_0 = vars(2);
delta = -vars(3);
epsilon = vars(4);

b = [E1_err^2,  0,       0,       0,      0,      0;
     0,       E2_err^2,  0,       0,      0,      0;
     0,       0,       E3_err^2,  0,      0,      0;
     0,       0,       0,       E4_err^2, 0,      0;
     0,       0,       0,       0,      E5_err^2, 0;
     0,       0,       0,       0,      0,      E6_err^2];
 
vars_errs = pinv(A)*b*pinv(A');
delta_1_err = sqrt(vars_errs(1,1));
delta_0_err = sqrt(vars_errs(2,2));
delta_err = sqrt(vars_errs(3,3));
epsilon_err = sqrt(vars_errs(4,4));

moment_ratio = abs(delta_0/delta_1);
moment_ratio_err = sqrt(delta_0_err^2*(1/delta_1)^2+...
                        delta_1_err^2*(delta_0/delta_1^2)^2);
H0 = delta_0/(2*mu_0);
H0_err = sqrt(delta_0_err^2*(1/(2*mu_0))^2+...
              mu_0_err^2*(delta_0/(2*mu_0^2))^2);
H1 = (3*delta_1)/(2*mu_1);
H1_err = sqrt(delta_1_err^2*(3/(2*mu_1))^2+...
              mu_1_err^2*((3*delta_1)/(2*mu_1^2))^2);

fprintf(['Delta_1 =   %.4e +/- %.1e eV\n'...
         'Delta_0 =   %.4e +/- %.1e eV\n'...
         'mu_0/mu_1 = %.4f +/- %.4f \n'...
         'delta   =   %.3e +/- %.1e eV\n'...
         'epsilon =   %.3e +/- %.1e eV\n'...
         'H_0 = %.2f +/- %.2f Tesla \n'...
         'H_1 = %.2f +/- %.2f Tesla \n'],...
         delta_1,delta_1_err,...
         delta_0, delta_0_err,...
         moment_ratio,moment_ratio_err,...
         delta, delta_err,...
         epsilon, epsilon_err,...
         H0, H0_err,...
         H1, H1_err)
     
s_delta_1 = get_sigma(delta_1, d_1_accepted, delta_1_err);
s_delta_0 = get_sigma(delta_0, d_0_accepted, delta_0_err);
s_moment_ratio = get_sigma(moment_ratio, moment_accepted, moment_ratio_err);
s_H_1 = get_sigma(H1, H_accepted, H1_err);
s_H_0 = get_sigma(H0, H_accepted, H0_err);

fprintf(['\n    sigma Delta_1 = %.2f \n'...
           '    sigma Delta_0 = %.2f \n'...
           '    sigma mu_0/mu_1 = %.2f \n'...
           '    sigma_H_0 = %.2f \n'...
           '    sigma_H_1 = %.2f \n\n'],...
            s_delta_1,...
            s_delta_0,...
            s_moment_ratio,...
            s_H_1,...
            s_H_0)

function s = get_sigma(obtained, accepted, obtained_err)
    p = normcdf(obtained, accepted, obtained_err);
    lower = p/2;
    upper = 1 - (p/2);
    bounds = norminv([lower upper]);
    s = bounds(2);
end







     
     