% Load in all the appropriate data
% Assumes main function has been ran to collected detected jumps
real_event_data = table2array(readtable("test_1xsnr_var_events.csv"));
real_times = real_event_data(:,1);
jumps_detected = jumps_detected_1xsnr_var;
jumps_measured = jumps_measured_1xsnr_var;
jump_stats = jump_stats_1xsnr_var;

% jumps_measured = table2array(readtable("10xSNR_measured.csv"));
% jumps_detected = table2array(readtable("10xSNR_detected.csv"));
% jump_stats = table2array(readtable("10xSNR_stats.csv"));

%% Jump Detection Evaluation
% Only evaluates our algorithms ability to detect any
% regardless if they are single or multi-events

labels = LabelDetectedJumps(jumps_detected, real_times, tmeas, tjump);
TP1 = length(jumps_detected(labels >= 1));
FP1 = length(jumps_detected(labels == 0));
FN1 = length(real_times) - TP1;
precision1 = TP1/(TP1+FP1);
recall1 = TP1/(TP1+FN1);

% beta = 2 will bias toward recall
% beta = 0.5 will bias precission
beta = 0.5;
Fscore1 =(1+beta)^2*(precision1*recall1)/(beta^2*precision1+recall1);


%% Jump Filtering Evaluation
% Evaluates our algorithms ability to detect any event after pre-filtering
% regardless if they are single or multi-events 
labels = LabelDetectedJumps(jumps_measured, real_times, tmeas, tjump);
TP2 = length(jumps_measured(labels >= 1));
FP2 = length(jumps_measured(labels == 0));
FN2 = length(real_times) - TP2;
precision2 = TP2/(TP2+FP2);
recall2 = TP2/(TP2+FN2);

% beta = 2 will bias toward recall
% beta = 0.5 will bias precission
beta = 0.5;
Fscore2 =(1+beta)^2*(precision2*recall2)/(beta^2*precision2+recall2);



%% Clustering/PostFiltering Process
% Note that Jump Stats and Jump Measures need to be defined at this point
% jumps_measured = jumps_measured_100xsnr_var;
% jump_stats = jump_stats_100xsnr_var;
select_features = [1,2,3,4,5];
desired_fraction = [.3];
eps_range = linspace(0.000,1,1000);

[final_clusters, final_epsilons, final_fracs]  = clustering(jump_stats, jumps_measured, desired_fraction, eps_range, select_features);
% Running clusterin gshould define the final_clusters variable
% used in next section
PostFilter = PostFiltering(jump_stats, 0);

%% Calculate Evaluation
% Note that at this point, real_times needs to be defined

% TP: detected single-jump event that are in the final cluster
% FP: Events in our final cluster that are either multi-event(label == 2)
% or no-event(label == 0)
for cluster_i = 0:2
    TP3 = length(jumps_measured(final_clusters >= cluster_i & labels == 1 & PostFilter == 1));
    FP3 = length(jumps_measured(final_clusters >= cluster_i & labels ~= 1 & PostFilter == 1));
    
    % FN: real-events that we did not detect as single-jumps
    FN3 = length(real_times) - TP3;
    precision3 = TP3/(TP3+FP3)
    recall3 = TP3/(TP3+FN3);
    
    % beta = 2 will bias toward recall
    % beta = 0.5 will bias precission
    beta = 0.5;
    Fscore =(1+beta)^2*(precision3*recall3)/(beta^2*precision3+recall3);
end