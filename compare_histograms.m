% Load in all the appropriate data
% Assumes main function has been with jump_stats, jumps_detected,
% jumps_measured stored in appropriate variables.
real_event_data = table2array(readtable("test_1xsnr_var_events.csv"));
real_times_1xsnr = real_event_data(:,1);
real_event_data = table2array(readtable("test_10xsnr_var_events.csv"));
real_times_10xsnr = real_event_data(:,1);
real_event_data = table2array(readtable("test_100xsnr_var_events.csv"));
real_times_100xsnr = real_event_data(:,1);

%% Detected Jump Labeling
% For each measured/detected event, label as 
% no-jump(0), single-jump(1) or multi-jump
labels = LabelDetectedJumps(jumps_measured_100xsnr_var, real_times_100xsnr, tmeas, tjump);

% choose two moments to plot labeling against
ind1=3;
ind2=8;
figure; hold on
% plot(jump_stats_1xsnr_var(:,ind1),jump_stats_1xsnr_var(:,ind2),'.','MarkerSize',7);
% plot(jump_stats_10xsnr_var(:,ind1),jump_stats_10xsnr_var(:,ind2),'.','MarkerSize',7);
plot(jump_stats_100xsnr_var(labels == 1,ind1),jump_stats_100xsnr_var(labels == 1,ind2),'.','MarkerSize',7);
plot(jump_stats_100xsnr_var(labels == 2,ind1),jump_stats_100xsnr_var(labels==2,ind2),'.','MarkerSize',7);
plot(jump_stats_100xsnr_var(labels == 0,ind1),jump_stats_100xsnr_var(labels==0,ind2),'.','MarkerSize',7);
plot(jump_stats_snrboost_fixed(:,ind1),jump_stats_snrboost_fixed(:,ind2),'.','MarkerSize',7);
legend('Variable event rate - single jumps','Variable event rate - multi jumps', 'Variable event rate - False Jumps','Fixed event rate - single jumps');
%% Histogram Plotting

figure;
for moment_i = 1:4
    subplot(2,2, moment_i)
    h1=histogram(jump_stats_100xsnr_var(labels==2,moment_i),100); hold on
    h2=histogram(jump_stats_100xsnr_var(labels==1,moment_i),100,'BinEdges',h1.BinEdges);
    legend('Variable event rate - multi jumps','Variable event rate - single jumps');
    pdist2(h1.BinCounts,h2.BinCounts)
end
% Plot of FWHM comparing multi-events to fixed rate-events
% figure;
% h1=histogram(jump_stats_100xsnr_var(labels==2,8),100); hold on
% h2=histogram(jump_stats_snrboost_fixed(:,8),100,'BinEdges',h1.BinEdges);
% legend('Variable event rate - multi jumps','Fixed event rate - single jumps');
% pdist2(h1.BinCounts,h2.BinCounts)




