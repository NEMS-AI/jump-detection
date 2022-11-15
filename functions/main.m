%% Load data
% Assume data is in correct format for csv
filenames = "isolated_1000snr.csv";
% filenames = "14_decay_exp_data.csv";

% Load in data and associated parameters
params = load_params(filenames);

% Iterate through time series data and calculate 
Fstats = calculate_fstats(params);
jumps_detected = detect(params, Fstats);
jumps_detected = preprocess(params, jumps_detected);

%TODO Move 
[jumps_measured, jump_stats, rel_jump_ts] = measure(params, jumps_detected, Fstats);


%% 

% TODO turn clustering into seperate function
% Figure what to do with current assortment of plots
% TODO Fix clustering for short datasets; uncomment to utilize
% clustering; 

select_features = [1,2];
desired_fraction = [.3];
eps_range = linspace(0.000,1,1000);

[final_clusters, final_epsilons, final_fracs]  = cluster_function(jump_stats, jumps_measured, desired_fraction, eps_range, select_features);