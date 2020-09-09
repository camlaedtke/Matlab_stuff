clear all
close all

data = dlmread('Time Measurements.csv',',',1,0);
t_nylon = data(:,1);
t_steel = data(:,2);
t_nylon_err = std(t_nylon)/5;
t_steel_err = std(t_steel)/5;

% k = (1:1:length(t_nylon))';
% fprintf('Row  t_nylon(s)  t_steel(s)     y(m)     y_err(m)\n');
% fprintf('%.0f      %.5f     %.5f     %.5f   %.5f\n', ...
%     [k,t_nylon,t_steel,y,y_err].')
data = dlmread('Ball Dropper Data.csv',',',2,5);
y = data(:,1)./1000; % drop heights (meters)
y_err = data(:,2)./1000; % meters
t_nylon = data(:,3);% seconds
t_steel = data(:,4);% seconds
m_nylon = 0.00409; % kg
m_steel = 0.02817; % kg

% delete = [4];
delete = [4,24];
% delete = [1,4,24];
% delete = [1,4,24,41];
% delete = [1,3,4,24,41];
[t_nylon,t_steel,y,y_err] = drop(delete,t_nylon,t_steel,y,y_err);

t_inf = (t_steel*m_steel- t_nylon*m_nylon)./(m_steel-m_nylon);
t_inf_err = 1/(m_steel-...
    m_nylon).*sqrt((m_steel.*t_steel_err).^2+....
    (m_nylon.*t_nylon_err).^2);


dydt_err = sqrt(t_inf_err^2.*(y./t_inf.^2).^2+y_err.^2.*(1./t_inf).^2);
fprintf(['\nTime error for nylon: %.0f us \nTime error for steel: %.0f us \n'...
 'Distance error: %.0f um \nTime error: %.0f us \n'...
 'y/t error: %.0f um/us to %.0f um/us\n'],t_nylon_err*10^6,t_steel_err*10^6,...
          y_err(1)*10^6,t_inf_err*10^6,min(dydt_err)*10^6,max(dydt_err*10^6));

del = (0.009:0.00002:0.014)';
sse = zeros(length(del),1);


for i = 1:length(del)
    del_0 = del(i);
    dydt = y./t_inf + del_0./t_inf;
    sse(i) = getSSE(t_inf,dydt,dydt_err,del_0);
   
end


figure(1)
set(gcf,'Units','Normalized','Position',[0.6 0.6 0.4 0.4])
scatter(del,sse,10)
xlabel('$\delta$ correction','Interpreter','latex')
ylabel('SSE','Interpreter','latex')
title('Correction factor vs LSQ SSE','Interpreter','latex')
set(gca,'FontSize',20) 
grid on

sse = -1.*sse;
[px,py] = findpeaks(sse,del);
del_best = py;
dydt = y./t_inf + del_best./t_inf;
p = getparams(t_inf,dydt,dydt_err,del_best);
yp = p.b.*t_inf + p.a;

figure(2)
set(gcf,'Units','Normalized','Position',[0.6 0.4 0.4 0.4])
scatter(t_inf,dydt,2)
hold on
plot(t_inf,yp)
xlabel('$t_{m \rightarrow \infty}$','Interpreter','latex')
ylabel('$\frac{y}{t_{m \rightarrow \infty}}$',...
    'Interpreter','latex')
title('Steel and Nylon Free Fall Fit (Correction)',...
    'Interpreter','latex')
set(gca,'FontSize',20) 
grid on
errorbar(t_inf,dydt,dydt_err,'o');
hold off


figure(3)
set(gcf,'Units','Normalized','Position',[0.6 0.0 0.4 0.4])
plot(t_inf(1:length(p.res),:), (p.res), 'x')
title('Weighted Residuals vs. x')
xlabel('x')
ylabel('Weighted Residuals')
set(gca,'FontSize',20) 
grid on


fprintf('\nBest Value of del is: del = %.5f\n',py)
fprintf('Sum of Chi Squared = %.0f \n',-px) 
fprintf('Reduced Chi Squared = %.1f \n',p.r_chi) 
fprintf('Calculated value of g = %.5f +/- %.5f m/s^2 \n',2*p.b,2*p.b_err)
fprintf('Calculated value of V0 = %.5f +/- %.5f m/s \n',p.a,p.a_err)
fprintf('\n ----------- FINISHED -------------\n');


function chival = getSSE(x,y,sig_y,del)
    [~, gof, ~] = fit(x, y,'poly1', 'Weights',(1./sig_y).^2);
    sprintf(['Chi Squared Fit --- del = %.4f'...
        ' --- SSE = %.0f\n'],del,gof.sse);
    chival = gof.sse;
end

function p = getparams(x,y,sig_y,del)
    [fitobj, gof, outp] = fit(x, y,'poly1', ...
        'Weights',(1./sig_y).^2);
    sprintf(['Chi Squared Fit --- del = %.5f'...
        ' --- SSE = %.0f\n'],del,gof.sse);
    p.b = fitobj.p1;
    p.a = fitobj.p2;
    p.res = outp.residuals;
    p.outp = outp;
    p.fitobj = fitobj;
    error_matrix = inv(outp.Jacobian'*outp.Jacobian); 
    dcov =  diag(error_matrix);
    unc = sqrt(dcov)';
    p.b_err = unc(1);
    p.a_err = unc(2);
    p.r_chi = gof.sse/gof.dfe;
end

function [t_nylon,t_steel,y,y_err] = drop(delete,t_nylon,t_steel,y,y_err)
    t_nylon(delete,:) = [];
    t_steel(delete,:) = [];
    y(delete,:) = [];
    y_err(delete,:) = [];
end


