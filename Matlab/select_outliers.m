%% Plot various jump statistics to find outliers

% alljumps = true(size(jumps_measured,1),1);
% pickjumps = false(size(jumps_measured,1),1);
% usejumps = alljumps;
% pickjumps = false(size(jumps_measured,1),1);
% % pickjumps = outliers;
% brush1 = [];
% brush2 = [];
% brush3 = [];
% brush4 = [];
% plotoutliers = 0;
usejumps = ~pickjumps;
% usejumps = pickjumps5;
% usejumps = alljumps;
plotoutliers = 0;

% usage 1:
% pickjumps1: manually selected high F stat max (points along fingerprint curve)
% pickjumps2: manually selected medium F stat max (points along fingerprint curve)
% pickjumps3: manually selected low F stat max (points along fingerprint curve)
% pickjumps4: manually selected very low F stat max (points along fingerprint curve)
% pickjumps5: manually selected outliers (points not along fingerprint curve)

% usage 2:
% usejumps: points along fingerprint curve (those not removed as outliers)
% pickjumps1: points with high standard deviation (sqrt var) categorized as outliers
% pickjumps2: points with low FWHM indicated as outliers

% usage 3:
% pickjumps6, pickjumps7: after selecting and plotting according "usage 1" or "usage 2", outliers
% are further sub-categorized according to additional selection criteria

% columns of jumps_measured:
% [Fstatmax t_mean_rel tvect(t_fwhm+1) tvect(t_bif+1) tvect(t_above_thresh+1) sqrt(F_var) peaks_right]

% co = [0.4660 0.6740 0.1880
%     0.8500 0.3250 0.0980
%     0 0.4470 0.7410];
% co = [0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840];
% set(groot,'defaultAxesColorOrder',co)

% scatter plot of fingerprint vector
figure; hold on
xlabel('Relative frequency shift (Mode 1)'); ylabel('Relative frequency shift (Mode 2)');
if plotoutliers == 1
    hold on
    plot(jumps_measured(pickjumps1,2),jumps_measured(pickjumps1,3),'.','MarkerSize',7);
    plot(jumps_measured(pickjumps2,2),jumps_measured(pickjumps2,3),'.','MarkerSize',7);
%     plot(jumps_measured(pickjumps3,2),jumps_measured(pickjumps3,3),'.','MarkerSize',7);
%     plot(jumps_measured(pickjumps4,2),jumps_measured(pickjumps4,3),'.','MarkerSize',7);
%     plot(jumps_measured(pickjumps5,2),jumps_measured(pickjumps5,3),'.','MarkerSize',7);
%     plot(jumps_measured(pickjumps6,2),jumps_measured(pickjumps6,3),'.','MarkerSize',7);
%     plot(jumps_measured(pickjumps7,2),jumps_measured(pickjumps7,3),'.','MarkerSize',7);
%     legend('Non-outliers (thus far)','High sqrt var','Low FWHM');
%     legend('Fingerprint curve','High sqrt var','Low FWHM','Remaining outliers');
%     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Large integral','Remaining outliers');
elseif plotoutliers == 2
    hold on
    plot(jumps_measured(pickjumps,2),jumps_measured(pickjumps,3),'.','MarkerSize',7);
end
plot(jumps_measured(usejumps,2),jumps_measured(usejumps,3),'.','MarkerSize',7);
% legend('High standard deviation','Low FWHM','Remaining data','Location','northwest');


% FWHM vs std dev
figure; hold on
% plot(jump_stats(usejumps,6),jump_stats(usejumps,3),'.','MarkerSize',7);
% set(gca, 'YScale', 'log');
if plotoutliers == 1
    hold on
    plot(jump_stats(pickjumps1,6),jump_stats(pickjumps1,3),'.','MarkerSize',7);
    plot(jump_stats(pickjumps2,6),jump_stats(pickjumps2,3),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps3,6),jump_stats(pickjumps3,3),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps4,6),jump_stats(pickjumps4,3),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps5,6),jump_stats(pickjumps5,3),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps6,6),jump_stats(pickjumps6,3),'.','MarkerSize',7);
%     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('Fingerprint curve','High Sqrt Var','Low FWHM','Remaining outliers');
end
plot(jump_stats(usejumps,6),jump_stats(usejumps,3),'.','MarkerSize',7);
xlabel('Standard deviation (s)');
ylabel('FWHM (s)');

%% Manually choose outliers or jumps of interest by interacting with figures
% Use figures to manually select jumps of interest, label them as "brush1", "brush2", etc.
% then flag corresponding indices among measured jumps and confirm selection

% % pick jumps on finger print scatter plot
% jumps_selected1 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush1,1)
%         if brush1(jj,1)==jumps_measured(ii,2) && brush1(jj,2)==jumps_measured(ii,3), jumps_selected1(ii)=1; end
%     end
% end
% % pick jumps on finger print scatter plot
% jumps_selected2 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush2,1)
%         if brush2(jj,1)==jumps_measured(ii,2) && brush2(jj,2)==jumps_measured(ii,3), jumps_selected2(ii)=1; end
%     end
% end
% % pick jumps on finger print scatter plot
% jumps_selected3 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush3,1)
%         if brush3(jj,1)==jumps_measured(ii,2) && brush3(jj,2)==jumps_measured(ii,3), jumps_selected3(ii)=1; end
%     end
% end
% % pick jumps on finger print scatter plot
% jumps_selected4 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush4,1)
%         if brush4(jj,1)==jumps_measured(ii,2) && brush4(jj,2)==jumps_measured(ii,3), jumps_selected4(ii)=1; end
%     end
% end
% % pick jumps on finger print scatter plot
% jumps_selected5 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush5,1)
%         if brush5(jj,1)==jumps_measured(ii,2) && brush5(jj,2)==jumps_measured(ii,3), jumps_selected5(ii)=1; end
%     end
% end

% pick jumps on Euclidean distance vs FWHM plot
% jumps_selected2 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush2,1)
%         if brush2(jj,1)==jump_stats(ii,5) && brush2(jj,2)==jump_stats_v_median(ii,1), jumps_selected2(ii)=1; end
%     end
% end

% jumps_selected6 = false(size(jumps_measured,1),1);
% for ii=1:length(usejumps)
%     for jj=1:size(brush1,1)
%         if brush1(jj,1)==jumps_measured(ii,2) && brush1(jj,2)==jumps_measured(ii,3), jumps_selected6(ii)=1; end
%     end
% end

% pickjumps = jumps_selected1 | jumps_selected2;

%% retrieve manually selected jumps
% jumps_selected are the manually selected points using the figure and "select data" features,
% these can be stored in the workspace and recovered into "pickjumps" using these sections.
% pickjumps1 = jumps_selected1;
% pickjumps2 = jumps_selected2;
% pickjumps3 = jumps_selected3;
% pickjumps4 = jumps_selected4;
% pickjumps5 = jumps_selected5;

%%
% pickjumps6 = jumps_selected6;

%% Plot median behavior ("jump signature") of various events of interest
% 
% [med_rel_jump0, med_Fstats0] = get_median_jump(alljumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% [med_rel_jump, med_Fstats] = get_median_jump(usejumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % plot_jump(nmodes,tvect(jump_range),med_rel_jump1,log10(med_Fstats1),log10(med_Fstats),tmeas,tjump,tjump_pre);
% % title('Large SNR');
% [med_rel_jump1, med_Fstats1] = get_median_jump(pickjumps1,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% [med_rel_jump2, med_Fstats2] = get_median_jump(pickjumps2,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % [med_rel_jump3, med_Fstats3] = get_median_jump(pickjumps3,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % [med_rel_jump4, med_Fstats4] = get_median_jump(pickjumps4,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % [med_rel_jump5, med_Fstats5] = get_median_jump(pickjumps5,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % [med_rel_jump6, med_Fstats6] = get_median_jump(pickjumps6,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% % [med_rel_jump7, med_Fstats7] = get_median_jump(pickjumps7,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% 
% % plot_jump(nmodes,tvect(jump_range),med_rel_jump2,log10(med_Fstats),log10(med_Fstats),tmeas,tjump,tjump_pre);
% % plot_jump(nmodes,tvect(jump_range),med_rel_jump1,log10(med_Fstats1),log10(med_Fstats0),tmeas,tjump,tjump_pre);
% % title('High Sqrt Var');
% % plot_jump(nmodes,tvect(jump_range),med_rel_jump2,log10(med_Fstats2),log10(med_Fstats0),tmeas,tjump,tjump_pre);
% % title('Low FWHM');
% plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),log10(med_Fstats0),tmeas,tjump,tjump_pre);
% % title('Remaining data');


[med_rel_jump, med_Fstats] = get_median_jump(pickjumps_mr,rel_jump_ts_1_pos,rel_jump_ts_2_pos,rel_jump_ts_3,Fstats_ts_pos);
plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),log10(med_Fstats),tmeas,tjump,tjump_pre);


%% Manually choose outliers or jumps of interest using cutoffs of statistics

% % pickjumps6 and 7 are optionally selected after a 1st round of determining outliers
% % to further classify remaining outliers
% pickjumps6 = false(size(jump_stats,1),1);
% pickjumps7 = pickjumps5;
% 
% for ji = 1:length(pickjumps6)
%     
%     Fint_left = jump_stats(ji,2);
%     if Fint_left > 20
%         pickjumps6(ji) = 1;
%         pickjumps7(ji) = 0;
%     end
% end

%% Manually choose outliers or jumps of interest using cutoffs of statistics

pickjumps0 = false(size(jump_stats,1),1);
pickjumps1 = false(size(jump_stats,1),1);
pickjumps2 = false(size(jump_stats,1),1);
pickjumps3 = false(size(jump_stats,1),1);
pickjumps4 = false(size(jump_stats,1),1);
pickjumps5 = false(size(jump_stats,1),1);
pickjumps = false(size(jump_stats,1),1);

for ji = 1:length(pickjumps)
    % columns of jump_stats:
    % [Fstatmax t_mean_rel tvect(t_fwhm+1) tvect(t_bif+1) tvect(t_above_thresh+1) sqrt(F_var) peaks_right]
    Fstatmax = jump_stats(ji,1);
    t_avg = jump_stats(ji,2);
    t_fwhm = jump_stats(ji,3);
    t_abovethresh = jump_stats(ji,5);
    t_sqrtvar = jump_stats(ji,6);
    rmse = jump_stats_v_median(ji,1);
    rmse_log = jump_stats_v_median(ji,2);
    pos_jump = jumps_measured(ji,2) > 0 || jumps_measured(ji,3) > 0;
    if nmodes ==3
        pos_jump = pos_jump || jumps_measured(ji,3) > 0;
    end
    
    % various versions of selection criteria commented out for reference, TODO: delete
%     if euc_dist > 30
%         pickjumps1(ji) = 1;
%     elseif t_fwhm < .06  % .0625 slightly better?
%         pickjumps2(ji) = 1;
%     elseif pos_jump
%         pickjumps3(ji) = 1;
%     elseif Fstatmax > 7E5   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
%         pickjumps4(ji) = 1;
%     elseif t_fwhm > .08
%         pickjumps5(ji) = 1;
%     end

%     if pos_jump
%         pickjumps0(ji) = 1;
%     elseif t_sqrtvar < .0125
%         pickjumps3(ji) = 1;
%     elseif rmse > .3
%         pickjumps1(ji) = 1;
%     elseif rmse > .135
%         pickjumps2(ji) = 1;
%     elseif t_sqrtvar > .018
%         pickjumps4(ji) = 1;
%     elseif Fstatmax > 7E5   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
%         pickjumps5(ji) = 1;
%     end

%     if pos_jump
%         pickjumps0(ji) = 1;
%     elseif rmse_log > .75
%         pickjumps1(ji) = 1;
%     elseif t_fwhm < .0365
%         pickjumps2(ji) = 1;
%     elseif Fstatmax > 7E5   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
%         pickjumps3(ji) = 1;
%     end
    
%     if pos_jump
%         pickjumps0(ji) = 1;
%     elseif t_sqrtvar > .034
%         pickjumps1(ji) = 1;
%     elseif t_fwhm < .03
%         pickjumps2(ji) = 1;
%     elseif Fstatmax > 1E6   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
%         pickjumps3(ji) = 1;
%     end
    
    if pos_jump
        pickjumps0(ji) = 1;
    elseif t_sqrtvar > .032
        pickjumps1(ji) = 1;
    elseif t_fwhm < .032
        pickjumps2(ji) = 1;
    elseif Fstatmax > 1E6   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
        pickjumps3(ji) = 1;
    end

end

pickjumps = pickjumps0 | pickjumps1 | pickjumps2; % do not include F stat max
% pickjumps = pickjumps0 | pickjumps1 | pickjumps2 | pickjumps3; % include F stat max
% pickjumps = pickjumps0; % include all negative jumps

%% Optionally plot specific jumps of interest. 

% % pickjumps = alljumps;
% % plotjumps = find(pickjumps);%[]; % list of indexes in jumps_measured to plot
% plotjumps = [16 25 42];
% 
% % optionally plot specific jumps of interest
% for pi = 1:length(plotjumps)
%     
%     ji = plotjumps(pi);
%     ti_jump = jumps_measured(ji,7); 
%     ti = ti_jump-Nmeas;
%     all_range = ti:ti+2*Nmeas+Njump-1;
%     rel_jump_ts = [rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)];
%     if nmodes==3
%         rel_jump_ts = [rel_jump_ts; rel_jump_ts_3(ji,:)];
%     end
%     plot_jump_simple(nmodes,tvect(all_range),rel_jump_ts,Fstats(all_range),tmeas,tjump,tjump_pre);
% end

%% pick out only positive jumps and match with cluster

% pickjumps_pos = ~pickjumps0;
% rel_jump_ts_1_pos = rel_jump_ts_1(pickjumps_pos,:);
% rel_jump_ts_2_pos = rel_jump_ts_2(pickjumps_pos,:);
% Fstats_ts_pos = Fstats_ts(pickjumps_pos,:);
% pickjumps_mr = logical(clusters(:,8));

