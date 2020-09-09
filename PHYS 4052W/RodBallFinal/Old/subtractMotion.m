function du_subtracted = subtractMotion(t, du)

    [fitobj, ~, ~] = fit(t', du', 'poly1');
    du_subtracted = du' - fitobj(t);
    
    
%     set(gcf,'units','normalized','position',[0.42 0.5 0.28 0.44],'MenuBar','none');
%     plot(t, du', 'LineWidth', 3, 'Color', 'r')
%     hold on
%     plot(t, fitobj(t),':', 'LineWidth', 4, 'Color', 'k')
%     grid on
%     title('Displacement vs. Time')
%     xlabel('Time (seconds)')
%     ylabel('Displacement (m)')


end

