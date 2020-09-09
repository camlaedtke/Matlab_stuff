function displacementresponse(t,du_right,v_right,t_data,du_data,v_data,params)
    figure(1)
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44],'MenuBar','none');
    plot(t, du_right, 'LineWidth', 3,'Color', 'r','DisplayName','right')
    hold on
    plot(t_data, du_data,':', 'LineWidth', 3, 'Color', 'k','DisplayName','Data')
    hold off
    title('Displacement vs. Time')
    xlabel('Time (seconds)','interpreter','latex')
    ylabel('Displacement (m)','interpreter','latex')
    xlim(params.tspan)
    grid on
    set(gca,'FontSize',params.font_size)
    legend

    figure(2)
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0.7 0.5 0.28 0.44],'MenuBar','none');
    patchline(t_data, v_data,'edgecolor','k','linewidth',params.linewidth,'edgealpha',params.edgealpha)
    hold on
    patchline(t, v_right,'edgecolor','r','linewidth',params.linewidth,'edgealpha',params.edgealpha)
    hold off
    title('Velocity vs. Time')
    xlabel('Time (seconds)','interpreter','latex')
    ylabel('Velocity (m/s)','interpreter','latex')
    xlim(params.tspan) 
    grid on
    set(gca,'FontSize',params.font_size)

end

