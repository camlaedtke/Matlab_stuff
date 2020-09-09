function displacementplots(t_sim, du_right, t_data, du_data,params)
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',params.fig_size);
    plot(t_sim, du_right, 'LineWidth', params.model_linewidth,'Color', 'r','DisplayName','Model')
    hold on
    plot(t_data, du_data,':', 'LineWidth', params.data_linewidth,'Color', 'k','DisplayName','Data')
    hold off
    title('Displacement at Free End')
    xlabel('Time [seconds]','interpreter','latex')
    ylabel('Displacement [meters]','interpreter','latex')
    xlim(params.t_range+params.xlim_shift) 
    grid on
    set(gca,'FontSize',params.font_size)
    legend
end

