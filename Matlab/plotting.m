% Mode 0 Reset Plotting
% Mode 1 Update outliers based on brushing 
% Mode 2 Use Thresholds for Outliers
% Mode 3 Display Selected Outlier Data
% Mode 4 Plot only Non-Outlier Data
% Mode 5 Write results to file

mode = 4;

if mode == 0
    %% Reset parameters used for outlier selection and display plots
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
    %% Update selection based on brush values on plot
    plotoutliers = 1;

    % Manual Outlier Detection     
    % pick jumps on finger print scatter plot
    jumps_selected1 = false(size(jumps_measured,1),1);
    for ii=1:length(usejumps)
        for jj=1:size(brush1,1)
            if brush1(jj,1)==jumps_measured(ii,2) && brush1(jj,2)==jumps_measured(ii,3), jumps_selected1(ii)=1; end
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
    %% Select outlisers based on rigidly defined thresholds
    plotoutliers = 1;
    ppickjumps = false(size(jump_stats,1),1);
    
    for ji = 1:length(pickjumps)
        
        Fstatmax = jump_stats(ji,1);
        t_abovethresh = jump_stats(ji,4);
        t_fwhm = jump_stats(ji,5);
        euc_dist = jump_stats_v_median(ji,1);
        pos_jump = jumps_measured(ji,2) > 0 || jumps_measured(ji,3) > 0;
        if nmodes ==3
            pos_jump = pos_jump || jumps_measured(ji,3) > 0;
        end
        
    %     if (euc_dist > 15 && Fstatmax > 1E5) || euc_dist > 100 || Fstatmax > 1E6 || t_fwhm < 0.06 || pos_jump  % 10x snr synthetic
        if euc_dist > 30 || Fstatmax > 1E4 || t_fwhm < 0.06 || pos_jump  % 1x snr synthetic
    %     if euc_dist > 30 || Fstatmax > 2E3 || t_fwhm < 0.06 || pos_jump  % 0.5x snr synthetic
          pickjumps(ji) = 1; 
        end
    end
    info_plots;
elseif mode == 3
    %% Optionally plot specific jumps of interest
    plotjumps = find(pickjumps);
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
    %% Reset parameters used for outlier selection and display plots
    plotoutliers = 0;

    [med_rel_jump, med_Fstats] = get_median_jump(~pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
    jump_range = 1:length(med_Fstats);
    plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),tmeas,tjump,tjump_pre);
    title('Non-outliers');
    info_plots;

elseif mode == 5
    %% save the usable jumps from jumps_measured (now having been processed)
    writematrix(jumps_measured(usejumps,:),uiputfile());
end

    
  


