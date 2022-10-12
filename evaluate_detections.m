% Load in all the appropriate data
% Assumes main function has been ran to collected detected jumps
real_event_data = table2array(readtable("test_10xsnr_var_events.csv"));
real_times = real_event_data(:,1);
% detected_times = jumps_measured(:,1);
jumps_measured = table2array(readtable("10xSNR_measured.csv"));
jumps_detected = table2array(readtable("10xSNR_detected.csv"));

jump_stats = table2array(readtable("10xSNR_stats.csv"));

%% Jump Detection Evaluation
detected_times = jumps_detected(:,1);
detected_starts = jumps_detected(:,2);
detected_ends = jumps_detected(:,3);
% Keep track of true positives and false positives
TP = 0;
FP = 0;
labeling = zeros(size(detected_times));

% Iterate through detected jumps and assign them as either a false positive
% or true negative
for detected_event_i = 1:length(detected_times)

    % Build a range around each detected jump
    detected_event = detected_times(detected_event_i);
    lower_range = detected_starts(detected_event_i);
    upper_range = detected_ends(detected_event_i);
    is_real = false;
    
    % Iterate through real events and see if any fall within jump range
    real_count = 0;
    for real_event_i=1:length(real_times)
        real_event = real_times(real_event_i);
        if lower_range <= real_event && real_event <= upper_range
            is_real = true;
            real_count = real_count + 1;
        end
    end

    % If matched with a real event, record positions and counts
    if is_real
        TP = TP + 1;
    else
        FP = FP + 1;
    end
    labeling(detected_event_i) = real_count;

end

FN = length(real_times) - sum(labeling);
[TP FP FN]


%% Jump Filtering Evaluation
detected_times = jumps_measured(:,1);
detected_starts = jumps_measured(:,2);
detected_ends = jumps_measured(:,3);

% Keep track of true positives and false positives
TP = 0;
FP = 0;
labeling = zeros(size(detected_times));

% Iterate through detected jumps and assign them as either a false positive
% or true negative
for detected_event_i = 1:length(detected_times)

    % Build a range around each detected jump
    detected_event = detected_times(detected_event_i);
    lower_range = detected_starts(detected_event_i);
    upper_range = detected_ends(detected_event_i);
    is_real = false;
    
    % Iterate through real events and see if any fall within jump range
    real_count = 0;
    for real_event_i=1:length(real_times)
        real_event = real_times(real_event_i);
        if lower_range <= real_event && real_event <= upper_range
            is_real = true;
            real_count = real_count + 1;
        end
    end

    % If matched with a real event, record positions and counts
    if is_real
        TP = TP + 1;
    else
        FP = FP + 1;
    end
    labeling(detected_event_i) = real_count;

end

FN = length(real_times) - sum(labeling);
[TP FP FN]



%% 
clustering;
PostFilter = PostFiltering(jump_stats, 75);

%%
single_event_count = length(detected_times(labeling == 1 & PostFilter == 1));
multi_event_count = length(detected_times(labeling > 1 & PostFilter == 1));
[single_event_count/(single_event_count+multi_event_count)]

cluster = 2;
for i = 1:cluster
    single_event_count_c = length(detected_times(final_clusters >= i & labeling == 1 & PostFilter == 1));
    multi_event_count_c = length(detected_times(final_clusters >= i & labeling > 1 & PostFilter == 1));
%     [single_event_count multi_event_count single_event_count_c multi_event_count_c]
    [single_event_count_c/(single_event_count_c+multi_event_count_c)]
end
