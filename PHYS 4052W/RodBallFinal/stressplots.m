function stressplots(t_sim, stress_sim, t_data, stress_data, params)
    % close all
    
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',params.fig_size);
    plot(t_sim, stress_sim,params.line_type,...
        'LineWidth', params.model_linewidth,'Color', 'r','DisplayName','Model')
    hold on
    plot(t_data, stress_data,':', 'LineWidth', params.data_linewidth,'Color', 'k','DisplayName','Data')
    hold off
    title('Stress at Free End')
    xlabel('Time [seconds]','interpreter','latex')
    ylabel('Stress [$N/m^{2}$]','interpreter','latex')
    xlim(params.t_range+params.xlim_shift) 
    ylim(params.ylim)
    grid on
    set(gca,'FontSize',params.font_size)
    legend
end

