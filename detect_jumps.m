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

% initialize vectors for selecting outliers
brush1 = [];
brush2 = [];
brush3 = [];
brush4 = [];

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
%     err_angle = 
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

alljumps = true(size(jumps_measured,1),1);
pickjumps = false(size(jumps_measured,1),1);

%% Calculate and plot median jump so far 

% usejumps = alljumps; % use all jumps (do this first)
usejumps = ~pickjumps; % use jumps except those selected as outliers according to criteria in next section

% all jumps
[med_rel_jump, med_Fstats] = get_median_jump(alljumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
jump_range = 1:length(med_Fstats);
plot_jump(nmodes,tvect(jump_range),med_rel_jump,log10(med_Fstats),log10(med_Fstats),tmeas,tjump,tjump_pre);
% title('All jumps');

% optionally plot specific jumps
% jump_range = 1:length(med_Fstats);
% for ji=1:size(jumps_selected6,1)
%     if jumps_selected6(ji)
%         jump_stats(ji,:)
%         plot_jump(nmodes,tvect(jump_range),[rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)],log10(Fstats_ts(ji,:)),log10(med_Fstats),tmeas,tjump,tjump_pre);
%     end
% end

%% Calculate deviation from median F stat vs time

tmax = tvect(med_Fstats==max(med_Fstats));
tthresh = tvect(med_Fstats==max(med_Fstats));
% use points with time in first 3/4 of pre-jump measurement window (up to
% "bifurcation point")
tweight1 = double(tvect(1:length(med_Fstats)) < .75*tmeas)';
% use points with time less than that of the maximum F stat
tweight2 = double(tvect(1:length(med_Fstats)) < tmax)';
% columns are (1) RMSE between F stat vs time and median such
% data, in the pre-jump measurement window, (2) RMSE of the log
% of that data.
jump_stats_v_median = zeros(size(jump_stats,1),2);
for ji = 1:size(jumps_measured,1)
    
    % various normalization methods applied to median statistics; some unused
    med_Fstats_norm = med_Fstats/max(med_Fstats);
    Fstats_ts_norm = Fstats_ts(ji,:)/max(Fstats_ts(ji,:));
    med_Fstats_log10 = log10(med_Fstats);
    Fstats_ts_log10 = log10(Fstats_ts(ji,:));
    med_Fstats_norm_log10 = log10(med_Fstats)/max(log10(med_Fstats));
    Fstats_ts_norm_log10 = log10(Fstats_ts(ji,:))/max(log10(Fstats_ts(ji,:)));
    rmse = sqrt(sum(tweight1.*(Fstats_ts_norm - med_Fstats_norm).^2)/sum(tweight1));
    rmse_log10 = sqrt(sum(tweight1.*(Fstats_ts_log10 - med_Fstats_log10).^2)/sum(tweight1));
    jump_stats_v_median(ji,1) = rmse;
    jump_stats_v_median(ji,2) = rmse_log10;
    
end

%% Plot various jump statistics to find outliers
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
% pickjumps2: manually selected medium F stat max (points along fingerprint curve)

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


% FWHM vs sqrt var
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
% legend('High standard deviation','Low FWHM','Remaining data');

% % t avg vs sqrt var
% figure;
% plot(jump_stats(usejumps,6),jump_stats(usejumps,2),'.','MarkerSize',7);
% xlabel('sqrt var (s)');
% ylabel('t avg (s)');
% % set(gca, 'YScale', 'log');
% if plotoutliers == 1
%     hold on
%     plot(jump_stats(pickjumps1,6),jump_stats(pickjumps1,2),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps2,6),jump_stats(pickjumps2,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps3,6),jump_stats(pickjumps3,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps4,6),jump_stats(pickjumps4,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps5,6),jump_stats(pickjumps5,2),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps6,6),jump_stats(pickjumps6,2),'.','MarkerSize',7);
% %     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('Fingerprint curve','High Sqrt Var','Low FWHM','Remaining outliers');
% end

% % rmse vs sqrt var
% figure;
% plot(jump_stats(usejumps,6),jump_stats_v_median(usejumps,2),'.','MarkerSize',7);
% xlabel('sqrt var (s)');
% ylabel('RMSE compared with median');
% % set(gca, 'YScale', 'log');
% if plotoutliers == 1
%     hold on
%     plot(jump_stats(pickjumps1,6),jump_stats_v_median(pickjumps1,2),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps2,6),jump_stats_v_median(pickjumps2,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps3,6),jump_stats_v_median(pickjumps3,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps4,6),jump_stats_v_median(pickjumps4,2),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps5,6),jump_stats_v_median(pickjumps5,2),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps6,6),jump_stats_v_median(pickjumps6,2),'.','MarkerSize',7);
% %     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('Fingerprint curve','High Sqrt Var','Low FWHM','Remaining outliers');
% end

% % FWHM vs F stat max
% figure;
% plot(jump_stats(usejumps,6),jump_stats(usejumps,1),'.','MarkerSize',7);
% xlabel('Sqrt var (s)');
% ylabel('F stat max');
% set(gca, 'YScale', 'log');
% if plotoutliers == 1
%     hold on
%     plot(jump_stats(pickjumps1,6),jump_stats(pickjumps1,1),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps2,6),jump_stats(pickjumps2,1),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps3,6),jump_stats(pickjumps3,3),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps4,6),jump_stats(pickjumps4,3),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps5,6),jump_stats(pickjumps5,3),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps6,6),jump_stats(pickjumps6,1),'.','MarkerSize',7);
% %     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('Fingerprint curve','High Sqrt Var','Low FWHM','Remaining outliers');
% end

% % Skew vs Kurt
% figure;
% plot(jump_stats(usejumps,6),jump_stats(usejumps,8),'.','MarkerSize',7);
% xlabel('Std dev');
% ylabel('Kurtosis');
% % set(gca, 'YScale', 'log');
% if plotoutliers == 1
%     hold on
%     plot(jump_stats(pickjumps1,6),jump_stats(pickjumps1,8),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps2,6),jump_stats(pickjumps2,8),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps3,7),jump_stats(pickjumps3,8),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps4,7),jump_stats(pickjumps4,8),'.','MarkerSize',7);
% %     plot(jump_stats(pickjumps5,7),jump_stats(pickjumps5,8),'.','MarkerSize',7);
%     plot(jump_stats(pickjumps6,6),jump_stats(pickjumps6,8),'.','MarkerSize',7);
% %     legend('High F stat max','Median F stat max','Low F stat max','Very Low F stat max','Outliers');
%     legend('Fingerprint curve','High Std dev','Low FWHM','Remaining outliers');
% end


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
pickjumps1 = jumps_selected1;
pickjumps2 = jumps_selected2;
pickjumps3 = jumps_selected3;
pickjumps4 = jumps_selected4;
pickjumps5 = jumps_selected5;

%%
pickjumps6 = jumps_selected6;

%% Plot median behavior of various events of interest
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


% below figures are used to plot median behavior of F stats vs time for
% various selected points of interest to help determine selection criteria
% for outliers

% figure;
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats1));%/max(log10(med_Fstats1)));
% hold on
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats2));%/max(log10(med_Fstats2)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats3));%/max(log10(med_Fstats3)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats4));%/max(log10(med_Fstats4)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats6));%/max(log10(med_Fstats6)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats7));%/max(log10(med_Fstats7)));
% 
% xlabel('Time (s)');
% ylabel('log(F statistic)');
% legend('High F stat max','Medium F stat max','Low F stat max','Very low F stat max','Large integral','Remaining outliers');
% 
% figure;
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats1)/max(log10(med_Fstats1)));
% hold on
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats2)/max(log10(med_Fstats2)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats3)/max(log10(med_Fstats3)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats4)/max(log10(med_Fstats4)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats6)/max(log10(med_Fstats6)));
% plot(tvect(jump_range)-tmeas-tjump_pre,log10(med_Fstats7)/max(log10(med_Fstats7)));
% 
% xlabel('Time (s)');
% ylabel('log(F statistic (normalized))');
% legend('High F stat max','Medium F stat max','Low F stat max','Very low F stat max','Large integral','Remaining outliers');
% 
% figure;
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats1/max(med_Fstats1));%/max(log10(med_Fstats1)));
% hold on
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats2/max(med_Fstats2));%/max(log10(med_Fstats2)));
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats3/max(med_Fstats3));%/max(log10(med_Fstats3)));
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats4/max(med_Fstats4));%/max(log10(med_Fstats4)));
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats6/max(med_Fstats6));%/max(log10(med_Fstats6)));
% plot(tvect(jump_range)-tmeas-tjump_pre,med_Fstats7/max(med_Fstats7));%/max(log10(med_Fstats7)));

% xlabel('Time (s)');
% ylabel('F statistic (normalized)');
% legend('High F stat max','Medium F stat max','Low F stat max','Very low F stat max','Large integral','Remaining outliers');
% % ylims = ylim;
% % ylim([ylims(1) ylims(2)*1.05]);
% % ylims = ylim;
% % xlims = xlim;
% % rectangle('Position',[-tjump_pre ylims(1) tjump ylims(2)-ylims(1)],'FaceColor', [1 0 0 0.1],'EdgeColor',[1 0 0 0]);
% % plot((-tmeas-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
% % plot((-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
% % plot((tjump-tjump_pre)*[1 1],ylim,'-','Color',[1 0 0 .8]);
% % plot((tjump-tjump_pre+tmeas)*[1 1],ylim,'-','Color',[1 0 0 .8]);

% [med_rel_jump2, med_Fstats2] = get_median_jump(pickjumps2,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% plot_jump(nmodes,tvect(jump_range),med_rel_jump2,log10(med_Fstats2),log10(med_Fstats),tmeas,tjump,tjump_pre);
% title('Medium SNR');

% [med_rel_jump3, med_Fstats3] = get_median_jump(pickjumps3,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% plot_jump(nmodes,tvect(jump_range),med_rel_jump3,log10(med_Fstats3),log10(med_Fstats),tmeas,tjump,tjump_pre);
% title('Small SNR');

% [med_rel_jump4, med_Fstats4] = get_median_jump(pickjumps4,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
% plot_jump(nmodes,tvect(jump_range),med_rel_jump4,log10(med_Fstats4),log10(med_Fstats),tmeas,tjump,tjump_pre);
% title('VSmall SNR');

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

%% Plot median behavior of various events of interest

[med_rel_jump1, med_Fstats1] = get_median_jump(pickjumps1,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump1,log10(med_Fstats1),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Large RMSE');

[med_rel_jump2, med_Fstats2] = get_median_jump(pickjumps2,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump2,log10(med_Fstats2),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Medium RMSE');

[med_rel_jump3, med_Fstats3] = get_median_jump(pickjumps3,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump3,log10(med_Fstats3),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Small variance');

[med_rel_jump4, med_Fstats4] = get_median_jump(pickjumps4,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump4,log10(med_Fstats4),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Large variance');

[med_rel_jump5, med_Fstats5] = get_median_jump(pickjumps5,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump5,log10(med_Fstats5),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Large F stat max');

[med_rel_jump0, med_Fstats0] = get_median_jump(~pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump0,log10(med_Fstats0),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Non-outliers');

%% Optionally plot specific jumps of interest. TODO: needs debugging

% pickjumps = alljumps;
% plotjumps = find(pickjumps);%[]; % list of indexes in jumps_measured to plot
plotjumps = [16 25 42];

% optionally plot specific jumps of interest
for pi = 1:length(plotjumps)
    
    ji = plotjumps(pi);
    ti_jump = jumps_measured(ji,7); 
    ti = ti_jump-Nmeas;
    all_range = ti:ti+2*Nmeas+Njump-1;
    rel_jump_ts = [rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)];
    if nmodes==3
        rel_jump_ts = [rel_jump_ts; rel_jump_ts_3(ji,:)];
    end
    plot_jump_simple(nmodes,tvect(all_range),rel_jump_ts,Fstats(all_range),tmeas,tjump,tjump_pre);
end

%% pick out only positive jumps and match with cluster

pickjumps_pos = ~pickjumps0;
rel_jump_ts_1_pos = rel_jump_ts_1(pickjumps_pos,:);
rel_jump_ts_2_pos = rel_jump_ts_2(pickjumps_pos,:);
Fstats_ts_pos = Fstats_ts(pickjumps_pos,:);
pickjumps_mr = logical(clusters(:,8));

%% Output chosen jumps
% have user choose output file to avoid overwriting
writematrix(jumps_measured(~pickjumps,:),uiputfile());
% writematrix(jumps_measured(usejumps,:),uiputfile());
writematrix(jump_stats(~pickjumps,:),uiputfile());
