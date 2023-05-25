function plot_jump_simple(nmodes,tvect,rel_jump_ts,Fstats,tmeas,tjump,tjump_pre)
% plot_jump plots relative jump for modes 1 and 2 and Fstats over
% measurement and jump intervals. tvect and Fstats inputs should be 
% passed in as the range of interest.
    showlegend = 1;
    normalize = 0;
    figure;
    plot(tvect-tmeas-tjump_pre,rel_jump_ts(1,:),'k'); hold on
    plot(tvect-tmeas-tjump_pre,rel_jump_ts(2,:),'b');
    if nmodes == 3
        plot(tvect-tmeas-tjump_pre,rel_jump_ts(3,:),'r');
    end
    xlabel('Time (s)');
    ylabel('Normalized jump');
    yyaxis right
    if normalize
        plot(tvect-tmeas-tjump_pre,Fstats/max(Fstats));
    else
        plot(tvect-tmeas-tjump_pre,Fstats);
    end
    ylabel('log(F statistic)');
%     ylabel('log(F statistic) (Normalized)');
    ylims = ylim;
    ylim([ylims(1) ylims(2)*1.05]);
    if showlegend
        if nmodes == 2
            legend('Mode 1','Mode 2','F statistic');
        elseif nmodes == 3
            legend('Mode 1','Mode 2','Mode 3','F statistic');
        end
    end
end

