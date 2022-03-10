function plot_jump(nmodes,tvect,rel_jump_ts,Fstats,Fstat_med,tmeas,tjump,tjump_pre)
% plot_jump plots relative jump for modes 1 and 2 and Fstats over
% measurement and jump intervals. tvect and Fstats inputs should be 
% passed in as the range of interest.
    showlegend = 0;
    figure;
    plot(tvect-tmeas-tjump_pre,rel_jump_ts(1,:),'k'); hold on
    plot(tvect-tmeas-tjump_pre,rel_jump_ts(2,:),'b');
    if nmodes == 3
        plot(tvect-tmeas-tjump_pre,rel_jump_ts(3,:),'r');
    end
    xlabel('Time (s)');
    ylabel('Normalized jump');
    yyaxis right
    plot(tvect-tmeas-tjump_pre,Fstats);
    plot(tvect-tmeas-tjump_pre,Fstat_med);
    ylabel('F statistic');
    ylims = ylim;
    xlims = xlim;
    rectangle('Position',[-tjump_pre ylims(1) tjump ylims(2)-ylims(1)],'FaceColor', [1 0 0 0.1],'EdgeColor',[1 0 0 0]);
    plot((-tmeas-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
    plot((-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
    plot((tjump-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
    plot((tjump-tjump_pre+tmeas)*[1 1],ylim,'-','Color',[1 0 0 .8]);
    if showlegend
        if nmodes == 2
            legend('Mode 1','Mode 2','F stat','F stat (median)','Jump window','Measurement windows');
        elseif nmodes == 3
            legend('Mode 1','Mode 2','Mode 3','F stat','F stat (median)','Jump window','Measurement windows');
        end
    end
end

