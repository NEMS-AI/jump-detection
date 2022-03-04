mode = 3;

if mode == 0
    % Reset parameters used for outlier selection and display plots
    brush1 = [];
    brush2 = [];
    brush3 = [];
    alljumps = true(size(jumps_measured,1),1);
    pickjumps = false(size(jumps_measured,1),1);

    usejumps = alljumps;
    plotoutliers = 0;

    [med_rel_jump, med_Fstats] = get_median_jump(~pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
    jump_range = 1:length(med_Fstats);
    plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),tmeas,tjump,tjump_pre);
    title('Non-outliers');
    info_plots;
    

elseif mode == 1
    plotoutliers = 1;

    % Manual Outlier Detection     
    % pick jumps on finger print scatter plot
    jumps_selected1 = false(size(jumps_measured,1),1);
    for ii=1:length(usejumps)
        for jj=1:size(brush1,1)
            if brush1(jj,1)==jumps_measured(ii,1) && brush1(jj,2)==jumps_measured(ii,2), jumps_selected1(ii)=1; end
        end
    end
    
    % pick jumps on F stat max vs time above threshold
    jumps_selected2 = false(size(jumps_measured,1),1);
    for ii=1:length(usejumps)
        for jj=1:size(brush2,1)
            if brush2(jj,1)==jump_stats(ii,4) && brush2(jj,2)==jump_stats(ii,1), jumps_selected2(ii)=1; end
        end
    end
    
    % pick jumps on FWHM vs time above threshold
    jumps_selected3 = false(size(jumps_measured,1),1);
    for ii=1:length(usejumps)
        for jj=1:size(brush3,1)
            if brush3(jj,1)==jump_stats(ii,4) && brush3(jj,2)==jump_stats(ii,5), jumps_selected3(ii)=1; end
        end
    end
    
    pickjumps = pickjumps|jumps_selected1 | jumps_selected2 | jumps_selected3;

    
    [med_rel_jump, med_Fstats] = get_median_jump(~pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
    jump_range = 1:length(med_Fstats);
    plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),tmeas,tjump,tjump_pre);
    title('Non-outliers');
    
    [med_rel_jump, med_Fstats] = get_median_jump(pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
    jump_range = 1:length(med_Fstats);
    plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),tmeas,tjump,tjump_pre);
    title('Outliers');
    info_plots;
elseif mode == 2
    plotoutliers = 1;
    pickjumps = false(size(jump_stats,1),1);
    for ji = 1:length(pickjumps)
        
        Fstatmax = jump_stats(ji,1);
        t_abovethresh = jump_stats(ji,4);
        t_fwhm = jump_stats(ji,5);
        npeaks = jump_stats(ji,end);
        
        if t_abovethresh > 0.2 || t_fwhm < 0.06 % 0.1x snr synthetic
    %     if t_abovethresh > 0.2 || t_abovethresh <0.14 || t_fwhm < 0.067 || npeaks > 1 % 1x snr
          pickjumps(ji) = 1; 
        end
    end
    info_plots;
elseif mode == 3
    %% Optionally plot specific jumps of interest
    
    plotjumps = find(pickjumps);%[]; % list of indexes in jumps_measured to plot
    
    % optionally plot specific jumps of interest
    for pi = 1:length(plotjumps)
        
        ji = plotjumps(pi);
        ti_jump = jumps_measured(ji,5); 
        ti = ti_jump-Nmeas;
        all_range = ti:ti+2*Nmeas+Njump-1;
        rel_jump_ts = [rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)];
        if nmodes==3
            rel_jump_ts = [rel_jump_ts; rel_jump_ts_3(ji,:)];
        end
        plot_jump(nmodes,tvect(all_range),rel_jump_ts,Fstats(all_range),tmeas,tjump,tjump_pre);
    end
elseif mode == 4
    % Should plot only the non-outlier data on with the info plots
    plotoutliers = 0;

end

    
  


