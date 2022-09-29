%% Set constants

nmodes = 2;
data_type = 'synthetic';        % choose from 'noise', 'data', 'synthetic'
jump_type = 'decay';            % choose from 'inst', 'decay'
if nmodes == 2
    tsample = .00025;           % seconds per sample
    tmeas = .2;          % time window for measuring jump (before and after samples)
    if strcmp(jump_type, 'inst')
        tjump = 0; 
        tjump_offset = 0;
    else
        % 10x SNR
        tjump = .06;            % time window for jump dynamic
        tjump_offset = -.005;     % jump time offset
        % 1x SNR
%         tjump = .1;            % time window for jump dynamic
%         tjump_offset = -.045;     % jump time offset
    end
elseif nmodes == 3
    tsample = 0.02;
    tmeas = 1;     % 1-10s optimal detection window
    tjump = 1;            % not sure exact PLL settings
end

% number of samples associated with different time windows
Nmeas = floor(tmeas/tsample); 
Njump = floor(tjump/tsample);
Njump_offset = floor(tjump_offset/tsample);

% F stat threshold for detection and measurement
% below corresponds to p value of 0.003 from bootstrapping for experimental dataset
Fstat_thresh = 490;

%% Load data

% TODO: MOVE PREPROCESSING CODE

[tvect, fvect] = load_data(filenames, data_type, nmodes);

% optionally select portion of data
if strcmp(data_type, 'data')
    if nmodes == 2
        tvectuse = tvect > 10 & tvect < 830;  % full dataset
        % tvectuse = tvect > 10 & tvect < 420;  % first half
        % tvectuse = tvect > 420 & tvect < 830;  % second half
        % tvectuse = tvect > 10 & tvect < 92;  % first tenth
        % tvectuse = tvect > 748 & tvect < 830;  % last tenth
        tvect = tvect(tvectuse);
        tvect = tvect-tvect(1);
        fvect = fvect(:,tvectuse);
        tstart = 10; % used for time lag calculation
        tlag_mode2_per_second = .01/830;
        tlag_mode3_per_second = 0;
        tlag_buffer = 100;
    elseif nmodes == 3
        tstart = 0; 
        tlag_mode2_per_second = 0;
        tlag_mode3_per_second = 0;
        tlag_buffer = 0;
    end
elseif strcmp(data_type, 'noise') || strcmp(data_type, 'synthetic')
    tstart = 0;
    tlag_mode2_per_second = 0;
    tlag_mode3_per_second = 0;
    tlag_buffer = 0;
    Nlag_m2 = 0;
    Nlag_m3 = 0;
end

% preprocess to remove drift? invert the below:
% f1_drift = (.99997-1)*mode1_fstart/700*tvect;
% f2_drift = (.99998-1)*mode2_fstart/500*tvect;

% last time point
ttot = length(tvect);
tfin = ttot-Nmeas*3-Njump;
% tfin = floor(tfin/30);   % optionally use subset of data to peek at results

% plot relative frequencies
figure;
plot(tvect,fvect(1,:)/fvect(1,1)-1,'k'); hold on
plot(tvect,fvect(2,:)/fvect(2,1)-1,'b');
if nmodes == 3
    plot(tvect,fvect(3,:)/fvect(3,1),'r');
end
xlabel('Time (s)')
ylabel('Relative frequency change');

%% Calculate F statistic vs time for entire dataset

Fstats = zeros(tfin,1);
fprintf('Fstats:         ');
for ti=1+tlag_buffer:tfin-tlag_buffer % allow for time lag in either direction
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 
    
    xi_range = ti:ti+Nmeas-1;
    yi_range = ti+Nmeas+Njump:ti+2*Nmeas+Njump-1;

    if tlag_mode2_per_second == 0 && (nmodes < 3 || tlag_mode3_per_second == 0)
        fvect_xi = fvect(:,xi_range);
        fvect_yi = fvect(:,yi_range);
    else
        get_freq_samples;
    end
    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    sigma_x = cov(fvect_xi');
    sigma_y = cov(fvect_yi');
    sigma_pool = sigma_x/2 + sigma_y/2;
    t2 = Nmeas/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
    p_dim = nmodes;
    Fstat = (2*Nmeas-p_dim-1)/(p_dim*(2*Nmeas-2))*t2;
    Fstats(ti+Nmeas) = Fstat;
    
end

% plot Fstats
figure;
plot(tvect(1:length(Fstats)),Fstats);
xlabel('Time (s)')
ylabel('F statistic');

%% Detect jumps

jumps_detected = [];
ti=1+tlag_buffer; % allow for time lag in either direction
fprintf('jumps_detected:         ');
while ti < tfin-tlag_buffer
    
    fprintf('\b\b\b\b\b\b\b%5.01f%%\n',ti/ttot*100); 
    if Fstats(ti) < Fstat_thresh
        ti = ti+1;
    else
        tii = 0;
        while Fstats(ti+tii) > Fstat_thresh && tii < tfin 
            tii = tii+1;
        end
        % columns are: Fstat max, time index of Fstat max, 
        % time index of peak estimated going backwards from right edge
        peakstats_simple = get_peak_stats_simple(Fstats,Fstat_thresh,ti,ti+tii); 
        Fstatmax = peakstats_simple(1);
        % by convention jump happens 1 sample before it can be measured
        % include manual offset in location of detected jump
        if strcmp(jump_type,'inst')
            ti_jump = peakstats_simple(2)+Njump_offset;
        else
            ti_jump = peakstats_simple(3)+Njump_offset;
        end
        % columns are time of detected jump, time of beginning of jump,
        % time of end of jump, time index of jump, ...
        % time index crossing above then below threshold, Fstat max
        jumps_detected = [jumps_detected; tvect(ti_jump) tvect(ti_jump)+tjump_offset ...
                          tvect(ti_jump)+tjump_offset+tjump ti_jump ti ti+tii Fstatmax ];
                      
        % start looking for next jump after this one
        ti = ti + tii + 1;
    end
end

%% Measure jumps that are tmeas apart and collect stats on those jumps

jumps_measured = [];
jump_stats = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];
rel_jump_ts_3 = [];
Fstats_ts = [];
fprintf('jumps_measured:         ');
% TODO: include first and last jumps
for ji = 2:size(jumps_detected,1)-1
    
    fprintf('\b\b\b\b\b\b\b%5.0f%%\n',ji/size(jumps_detected,1)*100); 
    t_prev = jumps_detected(ji-1,1);
    t_curr = jumps_detected(ji,1);
    t_next = jumps_detected(ji+1,1);
    t_begin = jumps_detected(ji,2); % start of jump (based on manual offset)
    t_end   = jumps_detected(ji,3); % end of jump (based on manual offset and tjump)
    ti_curr = jumps_detected(ji,4); % jump time
    t_above = jumps_detected(ji,5); % cross above threshold
    t_below = jumps_detected(ji,6); % cross below threshold
    Fstatmax = jumps_detected(ji,7);
    twidth_curr = tvect(t_below-t_above+1);

    if twidth_curr < tjump*2, continue; end
    
    % TODO: revisit this!
    if t_curr < t_prev + tmeas + 2*tjump || ...
       t_next < t_curr + tmeas + 2*tjump, continue; end

    ti = ti_curr-Nmeas;
    xi_range = ti:ti+Nmeas-1;
    yi_range = ti+Nmeas+Njump:ti+2*Nmeas+Njump-1;
    all_range = ti:ti+2*Nmeas+Njump-1;
    
    if tlag_mode2_per_second == 0 && (nmodes < 3 || tlag_mode3_per_second == 0)
        fvect_xi = fvect(:,xi_range);
        fvect_yi = fvect(:,yi_range);
    else
        get_freq_samples;
    end    

    xbar = median(fvect_xi,2);
    ybar = median(fvect_yi,2);
    sigma_x = cov(fvect_xi');
    sigma_y = cov(fvect_yi');
    sigma_pool = sigma_x/2 + sigma_y/2;
    rel_jump = ybar./xbar-1;
    rel_err = [sqrt(sigma_pool(1,1)) sqrt(sigma_pool(2,2))]./fvect(:,1)';
    err_rho = sigma_pool(1,2)/(sqrt(sigma_pool(1,1))*sqrt(sigma_pool(2,2)));

    % columns are: time of jump (s), time of beginning of jump, time of end of jump,
    % time index of jump, relative jumps (df1, df2, df3, etc.), relative error.
    jumps_measured = [jumps_measured; t_curr t_begin t_end ...
                                      ti_curr rel_jump' rel_err err_rho];

    peakstats = get_peak_stats(Fstats,tvect,t_above,t_below,tfin);  
    
    % columns are: Fstatmax, time average of F stats relative to jump time, F stat FWHM
    % standard deviation of F stat (for data above threshold)
    % skewness of F stat, kurtosis of F stat
    jump_stats = [jump_stats; peakstats];
    
    % Todo: use cells to improve flexbility
    if nmodes == 2
        rel_jump_ts = ([fvect(1,all_range); fvect(2,all_range+Nlag_m2)]-ybar)./(xbar-ybar);
    elseif nmodes == 3 
        rel_jump_ts = ([fvect(1,all_range); fvect(2,all_range+Nlag_m2); ...
                        fvect(3,all_range+Nlag_m3);]-ybar)./(xbar-ybar);
    end
    rel_jump_ts_1 = [rel_jump_ts_1; rel_jump_ts(1,:)];
    rel_jump_ts_2 = [rel_jump_ts_2; rel_jump_ts(2,:)];
    if nmodes == 3
        rel_jump_ts_3 = [rel_jump_ts_3; rel_jump_ts(3,:)];
    end
    Fstats_ts = [Fstats_ts; Fstats(all_range)'];
end

%% Output chosen jumps
% have user choose output file to avoid overwriting
writematrix(jumps_detected,uiputfile());
writematrix(jumps_measured,uiputfile());
writematrix(jump_stats,uiputfile());
writematrix(final_clusters,uiputfile());
