function PlotPMTData(params)

    clf
    fontsize = 20;
    f1 = get(groot, 'CurrentFigure');
    f1.Position = [100 100 800 800];
    e1 = errorbar(params.GmMax,params.xMax, params.xMax_err, 'o','DisplayName','Maxima');
    e1.Marker = 'o';
    e1.MarkerSize = 1;
    e1.Color = 'blue';
    e1.CapSize = 5;
    hold on
    plot(params.GmMax,params.Gm1,'b','DisplayName','Maxima Fit')
    e2 = errorbar(params.GmMin,params.xMin, params.xMin_err, 's','DisplayName','Minima');
    e2.Marker = 'o';
    e2.MarkerSize = 1;
    e2.Color = 'red';
    e2.CapSize = 5;
    plot(params.GmMin,params.Gm2,'r','DisplayName','Minima Fit')
    title('Plot of Peak Distances $X_{Maxima}$, $X_{Minima}$ vs. $G(m, k)$',...
        'FontSize',fontsize,'Interpreter','latex')
    xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
    ylabel('Relative Distance (meters)','FontSize',fontsize,'Interpreter','latex')
    legend('Location','northwest','FontSize',fontsize,'Interpreter','latex',...
        'NumColumns',2)
    hold off
    
    f2 = figure;
    hold on
    scatter(params.GmMax, params.max_chi,'filled');
    ylim([-3 3])
    f2.Position = [900 100 600 400];
    title('Weighted Residuas for Maxima Fit',...
        'FontSize',fontsize,'Interpreter','latex')
    xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
    grid on
    hold off
    
    f3 = figure;
    hold on
    scatter(params.GmMin,params.min_chi,'filled','r');
    ylim([-3 3])
    f3.Position = [900 500 600 400];
    title('Weighted Residuals for Minima Fit',...
        'FontSize',fontsize,'Interpreter','latex')
    xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
    grid on
    hold off
    
    hold on
    scatter(params.GmMin,params.min_chi,'filled','r');
    ylim([-3 3])
    f3.Position = [900 500 600 400];
    title('Weighted Residuals for Minima Fit',...
        'FontSize',fontsize,'Interpreter','latex')
    xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\chi$ value','FontSize',fontsize,'Interpreter','latex')
    grid on
    hold off
    
    f4 = figure;
    hold on
    scatter(params.GmMin,params.min_chi2,'filled','r');
    scatter(params.GmMax,params.max_chi2,'filled','b');
    f4.Position = [1500 100 800 600];
    title('Weighted $\chi^2$ values for Both Fits',...
        'FontSize',fontsize,'Interpreter','latex')
    xlabel('G(m)','FontSize',fontsize,'Interpreter','latex')
    ylabel('$\chi^2$ value','FontSize',fontsize,'Interpreter','latex')
    grid on
    hold off

end

