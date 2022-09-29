function [jumps_measured, jump_stats, rel_jump_ts_1] = measure(params, jumps_detected, Fstats)
% Measure jumps that are tmeas apart and collect stats on those jumps

% Variable Initialization
jumps_measured = [];
jump_stats = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];
Fstats_ts = [];
% TODO Clean up the use of params for initialization
tjump = params.tjump;
tmeas =  params.tmeas;
Fstat_thresh = params.Fstat_thresh;
Nmeas = params.Nmeas;
Njump = params.Njump;
fvect = params.fvect;
tfin = params.tfin;

% Iterate through each of the detected jumps 
fprintf('jumps_measured:         ');
% TODO: include first and last jumps
for ji = 2:size(jumps_detected,1)-1
    fprintf('\b\b\b\b\b\b\b%5.0f%%\n',ji/size(jumps_detected,1)*100);

    % Calculate time times of prev/current/next jumps
    t_prev = jumps_detected(ji-1,1);
    t_curr = jumps_detected(ji,1);
    ti_curr = jumps_detected(ji,3);
    t_next = jumps_detected(ji+1,1);
    Fstatmax = jumps_detected(ji,2);
    twidth_curr = params.tvect(jumps_detected(ji,5)-jumps_detected(ji,4)+1);

    % Skip detected jumps that are too short to likely be a real jump
    % and is a noise outlier instead
    % TODO: Check condition for validity and precision
    if twidth_curr < tjump, continue; end
    
    % Skip detected jumps that are too close to one another. The jump is
    % within another jump's measurement window.
    % TODO: Check condition for validity and precision
    if t_curr < t_prev + tmeas + tjump || ...
       t_next < t_curr + tmeas + tjump, continue; end


    ti_jump = jumps_detected(ji,3); 
    ti = ti_jump-Nmeas;
    
    xi_range = ti:ti+Nmeas-1;
    yi_range = ti+Nmeas+Njump:ti+2*Nmeas+Njump-1;
    all_range = ti:ti+2*Nmeas+Njump-1;
    
    fvect_xi = fvect(:,xi_range);
    fvect_yi = fvect(:,yi_range);

    % Calculate the frequency shift and associated error
    [rel_jump, rel_err, err_rho, xbar, ybar] = calculate_jump(params, fvect_xi, fvect_yi);


    % columns are: time of jump (s), df1, df2, df3, etc., relative error, 
    % time index of jump.
    jumps_measured = [jumps_measured; t_curr rel_jump' rel_err err_rho ti_curr];
    
    % TODO merge these times into get_peak_stats functions
    ti1 = xi_range(1);
    ti2 = xi_range(end);
    ti3 = all_range(end);
    ti4 = jumps_detected(ji,4);
    ti5 = jumps_detected(ji,5);
    peakstats = get_peak_stats(Fstats,params.tvect,ti1,ti2,ti3,ti4,ti5,tfin); 
    npeaks = count_peaks(Fstats,ti1,ti2,Fstat_thresh,Fstat_thresh,0.5,1.5); 
    
    % columns are: Fstatmax, time average of F stats relative to jump time, F stat FWHM
    % F stat time above bifurcation point (unused, set to 0), F stat time above threshold,
    % standard deviation of F stat (for data to the left of the jump and above threshold)
    % skewness of F stat, kurtosis of F stat, number of detected peaks (unused)
    jump_stats = [jump_stats; peakstats npeaks];
    
    % Get relative jump measurements for time series for each jump
    % Todo: use cells to improve flexbility
    jump_ts = [];
    for i = 1:params.nmodes
    jump_ts = [jump_ts; fvect(i,all_range)];
    end
    rel_jump_ts = (jump_ts-ybar)./(xbar-ybar);
    
end

    % Split up the relative frequencies for each mode
    rel_jump_ts_1 = [rel_jump_ts_1; rel_jump_ts(1,:)];
    rel_jump_ts_2 = [rel_jump_ts_2; rel_jump_ts(2,:)];

    Fstats_ts = [Fstats_ts; Fstats(all_range)'];
end
