%% LOOK FOR BAD POINTS
clear all
close all

data = dlmread('Time Measurements.csv',',',1,0);
t_nylon = data(:,1);
t_steel = data(:,2);
t_nylon_err = std(t_nylon);
t_steel_err = std(t_steel);


% plot steel and nylon by themselves to spot outliers
% Different Heights
data = dlmread('Ball Dropper Data.csv',',',2,5);
y = data(:,1)./1000; % drop heights (meters)
y_err = data(:,2)./1000; % meters
t_nylon = data(:,3);% seconds
t_steel = data(:,4);% seconds

m_nylon = 0.00409; % kg
m_steel = 0.02817; % kg

t_inf = (t_steel*m_steel- t_nylon*m_nylon)./(m_steel-m_nylon);
t_inf_err = (1./(m_steel-m_nylon)).*sqrt((m_steel.*t_steel_err).^2+....
    (m_nylon.*t_nylon_err).^2);
t_inf_err = t_inf_err.*ones(length(t_inf),1);

dydt = y./t_inf; 
dydt_err = sqrt(y_err.^2.*(1./t_inf).^2 + t_inf_err.^2.*(y./t_inf.^2).^2);

analyzedata(t_inf,t_inf_err,dydt_err,y,y_err)

% delete = [4,24];
% [t_inf,t_inf_err,dydt,dydt_err,y,y_err] = droppoints(delete,t_inf,t_inf_err,dydt,dydt_err,y,y_err);
% analyzedata(t_inf,t_inf_err,dydt_err,y,y_err)

% delete = [1,41];
% [t_inf,dydt,dydt_err,y,y_err] = droppoints(delete,t_inf, dydt,dydt_err,y,y_err);
% analyzedata(t_inf,dydt_err,y)

% delete = [1,4,41];
% [t_inf,dydt,dydt_err,y,y_err] = droppoints(delete,t_inf, dydt,dydt_err,y,y_err);
% analyzedata(t_inf,dydt_err,y)
% 
% delete = [1,3,4,41];
% [t_inf,dydt,dydt_err,y,y_err] = droppoints(delete,t_inf, dydt,dydt_err,y,y_err);
% analyzedata(t_inf,dydt_err,y)
% 
% delete = [4,19,20,22,23,24];
% [t_inf,t_inf_err,dydt,dydt_err,y,y_err] = droppoints(delete,t_inf,t_inf_err, dydt,dydt_err,y,y_err);
% analyzedata(t_inf,t_inf_err,dydt_err,y,y_err)


%% FUNCTION

function [t_inf,t_inf_err, dydt, dydt_err, y, y_err] = droppoints(delete,t_inf,t_inf_err,dydt,dydt_err,y,y_err)
    t_inf(delete,:) = [];
    t_inf_err(delete,:) = [];
    dydt(delete,:) = [];
    dydt_err(delete,:) = [];
    y(delete,:) = [];
    y_err(delete,:) = [];
end

function chival = getSSE(x,y,sig_y,del)
    [~, gof, ~] = fit(x, y,'poly1', 'Weights',(1./sig_y).^2);
    sprintf('Chi Squared Fit --- del = %.4f --- SSE = %.3f\n',del,gof.sse);
    chival = gof.sse;
end

function p = getparams(x,y,sig_y,del)
    [fitobj, gof, outp] = fit(x, y,'poly1', 'Weights',(1./sig_y).^2);
    fprintf('Chi Squared Fit --- del = %.5f --- SSE = %.3f\n',del,gof.sse);
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
end

function analyzedata(t_inf,t_inf_err,dydt_err,y,y_err)
    del = (0.006:0.00001:0.018)';
    sse = zeros(length(del),1);
      
    for i = 1:length(del)
        del_0 = del(i);
        dydt = y./t_inf + del_0./t_inf;
        dydt_err = sqrt(y_err.^2.*(1./t_inf).^2 + t_inf_err.^2.*(y./t_inf.^2).^2);
        sse(i) = getSSE(t_inf,dydt,dydt_err,del_0);
    end
    
    figure(1)
    set(gcf,'Units','Normalized','Position',[0.05 0.5 0.3 0.4])
    scatter(del,sse,10)
    xlabel('b'' correction','Interpreter','latex')
    ylabel('SSE','Interpreter','latex')
    title('Correction factor vs LSQ SSE','Interpreter','latex')
    set(gca,'FontSize',20) 
    grid on
    
    sse = -1.*sse;
    [px,py] = findpeaks(sse,del);
    fprintf('\n\nBest Value of del is: del = %.5f\n',py)
    fprintf('Minimum Sum of Chi Squared + 1 = %.5f \n',-px+1) 
    del_best = py;
    dydt = y./t_inf + del_best./t_inf;
  
    p = getparams(t_inf,dydt,dydt_err,del_best);
    yp = p.b.*t_inf + p.a;
    fprintf('Calculated value of g = %.5f +/- %.5f m/s^2 \n',2*p.b,2*p.b_err)
    fprintf('Calculated value of V0 = %.5f +/- %.5f m/s \n\n',p.a,p.a_err)
    
    figure(2)
    set(gcf,'Units','Normalized','Position',[0.35 0.5 0.4 0.5])
    scatter(t_inf,dydt,5)
    hold on
    plot(t_inf,yp)
    xlabel('$t_{m \rightarrow \infty}$','Interpreter','latex')
    ylabel('$\frac{y}{t_{m \rightarrow \infty}}$','Interpreter','latex')
    title('Steel and Nylon Free Fall Fit (Correction)','Interpreter','latex')
    set(gca,'FontSize',20) 
    grid on
    errorbar(t_inf,dydt,dydt_err,'o');
    hold off
    
    % Look for bad points
    figure(3)
    set(gcf,'Units','Normalized','Position',[0.65 0.5 0.5 0.5])
    plot(t_inf(1:length(p.res),:), (p.res), 'x')
    title('Weighted Residuals vs. x')
    xlabel('x')
    ylabel('Weighted Residuals')
    set(gca,'FontSize',20) 
    grid on
    
   
end






