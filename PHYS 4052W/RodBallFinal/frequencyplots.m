function frequencyplots(data_freqs, data_FFT, sim_freqs, sim_FFT,params)
    %close all
    %figure(1)
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',params.fig_size,'MenuBar','none');
    plot(data_freqs, data_FFT,'Color','k','linewidth',params.data_linewidth)
    hold on
    plot(sim_freqs, sim_FFT,'Color','r','linewidth',params.model_linewidth)
    hold off
    title(params.title)
    xlabel(params.xlabel,'interpreter','latex')
    ylabel('Amplitude','interpreter','latex')
    xlim(params.freq_range + params.freq_shift)
    ylim(params.ylim)
    grid on
    legend('Model', 'Data')
    set(gca,'FontSize',params.font_size,'YScale', params.yscale)
end

