function D = odebifurcationplot(params, conds)
    clf
    conds.numpoints = conds.tend - conds.bifurcationstart;
    % range of gamma values to sample
    gamma = linspace(conds.gstart,conds.gend,...
             (conds.gend-conds.gstart)*conds.samplefreq+1);
    rangegamma = length(gamma);
    D.x = zeros(conds.numpoints, rangegamma);
    D.y = zeros(conds.numpoints, rangegamma);

    p = scatter(D.x(:), D.y(:),conds.pointsize,...
        'XDataSource','D.x(:)','YDataSource', 'D.y(:)');

    fig = get(groot,'CurrentFigure');
    fig.Position = [200 200 1000 500];
    ax1 = gca;
    ax1.Position = [0.1 0.1 0.8 0.8];
    ax1.Color = 'White';
    ax1.LabelFontSizeMultiplier = 1.5;
    ax1.TitleFontSizeMultiplier = 1.75;
    ax1.XLim = [conds.gstart conds.gend];
    ax1.YLim = [-1 0];
    ax1.XAxisLocation = 'bottom';
    ax1.XLabel.String = 'Gamma Value';
    ax1.YLabel.String = 'Theta';
    title('Bifurcation Diagram')
    xlabel('gamma')
    ylabel('theta')
    linkdata 'on'

    pb = CmdLineProgressBar('Progress ...');
    for i = 1:rangegamma
        params.gamma = gamma(i);
        conds.tstart = 0;
        equ = odegetddpequation(params);
        sol = odenumericsolve(equ, conds);
        conds.tstart = conds.bifurcationstart;
        S = odepoincaredata(sol, conds);

        D.x(:,i) = gamma(i);
        D.y(:,i) = S.x;   

        refreshdata
        drawnow

        pb.print(i,rangegamma)
    end
end

