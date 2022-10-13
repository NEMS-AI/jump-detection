% Load in all the appropriate data
% Assumes main function has been with jump_stats, jumps_detected,
% jumps_measured stored in appropriate variables.
real_event_data = table2array(readtable("test_1xsnr_var_events.csv"));
real_times_1xsnr = real_event_data(:,1);
real_event_data = table2array(readtable("test_10xsnr_var_events.csv"));
real_times_10xsnr = real_event_data(:,1);
real_event_data = table2array(readtable("test_100xsnr_var_events.csv"));
real_times_100xsnr = real_event_data(:,1);

%% Jump Detection Evaluation
% Keep track of double jumps
single_jump_100xsnr_var = false(size(jumps_measured_100xsnr_var(:,1)));
multi_jump_100xsnr_var = false(size(single_jump_100xsnr_var));
fp_jump_100xsnr_var = false(size(single_jump_100xsnr_var));

% Iterate through detected jumps and assign them as either a false positive
% or true negative
for measured_event_i = 1:size(jumps_measured_100xsnr_var,1)

    % Build a range around each detected jump
    detected_event = jumps_measured_100xsnr_var(measured_event_i,1);
    lower_range = jumps_measured_100xsnr_var(measured_event_i,2);
    upper_range = jumps_measured_100xsnr_var(measured_event_i,3);
    
    % Iterate through real events and see if any fall within jump range
    for real_event_i=1:length(real_times_100xsnr)
        real_event = real_times_100xsnr(real_event_i);
        if lower_range <= real_event && real_event <= upper_range
            single_jump_100xsnr_var(measured_event_i) = true;
        end
    end

    real_count = 0;
    for real_event_i=1:length(real_times_100xsnr)
        real_event = real_times_100xsnr(real_event_i);
        if lower_range-tmeas-tjump <= real_event && real_event <= upper_range+tmeas+tjump
            real_count = real_count + 1;
            if real_count > 1
                multi_jump_100xsnr_var(measured_event_i) = true;
                single_jump_100xsnr_var(measured_event_i) = false;
            end
        end
    end
end

ind1=3;
ind2=8;
figure; hold on
% plot(jump_stats_1xsnr_var(:,ind1),jump_stats_1xsnr_var(:,ind2),'.','MarkerSize',7);
% plot(jump_stats_10xsnr_var(:,ind1),jump_stats_10xsnr_var(:,ind2),'.','MarkerSize',7);
plot(jump_stats_100xsnr_var(single_jump_100xsnr_var,ind1),jump_stats_100xsnr_var(single_jump_100xsnr_var,ind2),'.','MarkerSize',7);
plot(jump_stats_100xsnr_var(multi_jump_100xsnr_var,ind1),jump_stats_100xsnr_var(multi_jump_100xsnr_var,ind2),'.','MarkerSize',7);
plot(jump_stats_snrboost_fixed(:,ind1),jump_stats_snrboost_fixed(:,ind2),'.','MarkerSize',7);
legend('Variable event rate - single jumps','Variable event rate - multi jumps','Fixed event rate - single jumps');

figure;
h1=histogram(jump_stats_100xsnr_var(multi_jump_100xsnr_var,ind1),100); hold on
h2=histogram(jump_stats_snrboost_fixed(:,ind1),100,'BinEdges',h1.BinEdges);
legend('Variable event rate - multi jumps','Fixed event rate - single jumps');
pdist2(h1.BinCounts,h2.BinCounts)