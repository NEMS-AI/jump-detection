function params = load_params(filenames)
%% Set constants
params.jump_type = "decay";
params.nmodes = 2;
params.data_type = 'synthetic';     % choose from 'noise', 'data', 'synthetic'
params.tsample = .00025;        % seconds per sample
params.tmeas_detect = .1;      % time window for detecting jump (before and after samples)
params.tjump = .09;             % time window for (full) jump itself; for 1-10x snr
params.tjump_pre = .08;         % portion of jump window before F threshold crossing; for 1-10x snr
%     tjump = .12;             % time window for (full) jump itself; for 0.5x snr
%     tjump_pre = .11;         % portion of jump window before F threshold crossing; for 0.5x snr
params.tmeas = .2;              % 200 ms for final measurement of jump height

params.tjump_post = params.tjump - params.tjump_pre; % portion of jump window after F threshold crossing

% number of samples associated with different time windows
params.Nmeas = floor(params.tmeas/params.tsample); 
params.Njump = floor(params.tjump/params.tsample);

% F stat threshold for detection and measurement
params.Fstat_thresh = 700;

%% Read file
syn_data = readtable(filenames(1));
syn_data = table2array(syn_data);
tvect = syn_data(:,1);
params.tvect=tvect-tvect(1);
fvect = syn_data(:,2:2+params.nmodes-1)';
params.fvect = fvect;

params.tstart = 0;
params.tlag_mode2_per_second = 0;
params.tlag_mode3_per_second = 0;
params.tlag_buffer = 0;
params.Nlag_m2 = 0;
params.Nlag_m3 = 0;

params.ttot = length(params.tvect);
params.tfin = params.ttot-params.Nmeas*3-params.Njump;


params.tjump = .06;            % time window for jump dynamic
params.tjump_offset = -.005;     % jump time offset
params.Njump_offset = floor(params.tjump_offset/params.tsample);