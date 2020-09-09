function frequencyplot(freqs, FFT, params)
    close all
    figure(1)
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0.42 0.5 0.56 0.44])%,'MenuBar','none');
    plot(freqs, FFT,'Color','k','linewidth',params.data_linewidth)
    title('Fourier Transform of the Data and Model (Velocity)')
    xlabel(params.xlabel,'interpreter','latex')
    ylabel('Amplitude','interpreter','latex')
    xlim(params.freq_range + params.freq_shift)
    ylim(params.ylim)
    grid on
    legend('Model', 'Data')
    set(gca,'FontSize',params.font_size,'YScale', params.yscale)
end

