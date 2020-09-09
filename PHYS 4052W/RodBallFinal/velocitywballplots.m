function velocitywballplots(t_sim, v_left, v_ball, params)
    close all
    figure(2)
    set(gcf,'units','normalized','position',[0.42 0.04 0.56 0.44],'MenuBar','none');
    plot(t_sim*1000, v_left, 'LineWidth', params.model_linewidth,'Color', 'b','DisplayName','First Mass')
    hold on
    plot(t_sim*1000, v_ball,'LineWidth', params.data_linewidth,'Color', 'g','DisplayName','Ball')
    hold off
    title('Model: Velocity of Ball and First Mass at Contact End')
    xlabel('Time [msec]','interpreter','latex')
    ylabel('Velocity [m/s]','interpreter','latex')
    xlim(params.t_range+params.xlim_shift) 
    grid on
    set(gca,'FontSize',params.font_size)
    legend
end

