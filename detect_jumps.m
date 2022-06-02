%% Set constants

nmodes = 2;
data_type = 'synthetic';     % choose from 'noise', 'data', 'synthetic'
if nmodes == 2
    tsample = .00025;        % seconds per sample
    tmeas_detect = .1;      % time window for detecting jump (before and after samples)
    tjump = .09;             % time window for (full) jump itself; for 1-10x snr
    tjump_pre = .08;         % portion of jump window before F threshold crossing; for 1-10x snr
%     tjump = .12;             % time window for (full) jump itself; for 0.5x snr
%     tjump_pre = .11;         % portion of jump window before F threshold crossing; for 0.5x snr
    tmeas = .2;              % 200 ms for final measurement of jump height
elseif nmodes == 3
    tsample = 0.02;
    tmeas_detect = 1;     % 1-10s optimal detection window
    tjump = 1;            % not sure exact PLL settings
    tjump_pre = 0.5;
    tmeas = 1; 
end
tjump_post = tjump - tjump_pre; % portion of jump window after F threshold crossing

% number of samples associated with different time windows
Ndetect = floor(tmeas_detect/tsample); 
Njump = floor(tjump/tsample);
Npre = floor(tjump_pre/tsample);

% F stat threshold for detection and measurement
Fstat_thresh_detect = 300;
Fstat_thresh_meas = 600;

%% Load data

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
tfin = ttot-Ndetect*3-Njump;
% tfin = floor(tfin/30);   % optionally use subset of data to peek at results

% plot relative frequencies
figure;
plot(tvect,fvect(1,:)/fvect(1,1),'k'); hold on
plot(tvect,fvect(2,:)/fvect(2,1),'b');
if nmodes == 3
    plot(tvect,fvect(3,:)/fvect(3,1),'r');
end
xlabel('Time (s)')
ylabel('Relative frequency');

%% Calculate F statistic vs time for entire dataset

Fstats = zeros(tfin,1);
fprintf('Fstats:         ');
for ti=1+tlag_buffer:tfin-tlag_buffer % allow for time lag in either direction
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 
    
    xi_range = ti:ti+Ndetect-1;
    yi_range = ti+Ndetect+Njump:ti+2*Ndetect+Njump-1;

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
    t2 = Ndetect/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
    p_dim = nmodes;
    Fstat = (2*Ndetect-p_dim-1)/(p_dim*(2*Ndetect-2))*t2;
    Fstats(ti+Ndetect) = Fstat;
    
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
    if Fstats(ti) < Fstat_thresh_detect
        ti = ti+1;
    else
        tii = 0;
        while Fstats(ti+tii) > Fstat_thresh_detect && tii < tfin 
            tii = tii+1;
        end
        t_above_thresh = tii;
        Fstatmax = max(Fstats(ti:ti+tii));
        % by convention jump happens 1 sample before it can be measured
        ti_jump = ti+tii-Npre-1;
        % columns are time of detected jump, max of F statistic of jump, 
        % time index of jump, time index crossing above then below threshold
        jumps_detected = [jumps_detected; tvect(ti_jump) Fstatmax ti_jump ti ti+tii];
                      
        % start looking for next jump after this one
        ti = ti + tii + 1;
    end
end

%% Measure jumps that are tmeas apart and collect stats on those jumps

Nmeas = floor(tmeas/tsample);
jumps_measured = [];
jump_stats = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];
rel_jump_ts_3 = [];
Fstats_ts = [];

fprintf('jumps_measured:         ');
for ji = 2:size(jumps_detected,1)-1
        
    fprintf('\b\b\b\b\b\b\b%5.0f%%\n',ji/size(jumps_detected,1)*100); 
    t_prev = jumps_detected(ji-1,1);
    t_curr = jumps_detected(ji,1);
    ti_curr = jumps_detected(ji,3);
    t_next = jumps_detected(ji+1,1);
    Fstatmax = jumps_detected(ji,2);
    
    if t_curr - tjump_pre < t_prev + tmeas + tjump_post || ...
       t_next - tjump_pre < t_curr + tmeas + tjump_post, continue; end

    if Fstatmax < Fstat_thresh_meas, continue; end

    ti_jump = jumps_detected(ji,3); 
    ti = ti_jump-Nmeas;
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

    % columns are: time of jump (s), df1, df2, df3, etc., relative error, 
    % time index of jump.
    jumps_measured = [jumps_measured; t_curr rel_jump' rel_err err_rho ti_curr];
    
    ti1 = xi_range(1);
    ti2 = xi_range(end);
    ti3 = all_range(end);
    ti4 = jumps_detected(ji,4);
    ti5 = jumps_detected(ji,5);
    peakstats = get_peak_stats(Fstats,tvect,ti1,ti2,ti3,ti4,ti5,tfin); 
    npeaks = count_peaks(Fstats,ti1,ti2,Fstat_thresh_detect,Fstat_thresh_meas,0.5,1.5); 
    
    % columns are: Fstatmax, time average of F stats relative to jump time, F stat FWHM
    % F stat time above bifurcation point (unused, set to 0), F stat time above threshold,
    % standard deviation of F stat (for data to the left of the jump and above threshold)
    % skewness of F stat, kurtosis of F stat, number of detected peaks (unused)
    jump_stats = [jump_stats; peakstats npeaks];
    
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
writematrix(jumps_measured,uiputfile());
writematrix(jump_stats,uiputfile());
