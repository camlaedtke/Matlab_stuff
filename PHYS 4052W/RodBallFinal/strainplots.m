function strainplots(t_sim, strain_sim_right, t_data, strain_data, params)
    % close all
    
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',params.fig_size);
    plot(t_sim, strain_sim_right,params.line_type,...
        'LineWidth', params.model_linewidth,'Color', 'r','DisplayName','Model')
    hold on
    plot(t_data, strain_data,':', 'LineWidth', params.data_linewidth,'Color', 'k','DisplayName','Data')
    hold off
    title('Strain at Free End')
    xlabel('Time [seconds]','interpreter','latex')
    ylabel('Strain','interpreter','latex')
    ylim([params.ylim])
    xlim(params.t_range+params.xlim_shift) 
    grid on
    set(gca,'FontSize',params.font_size)
    legend
end

