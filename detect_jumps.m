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
%     err_angle = 
    % columns are: time of jump (s), df1, df2, df3, etc.
    jumps_measured = [jumps_measured; t_curr rel_jump' rel_err];
    
    ti1 = jumps_detected(ji,4);
    ti2 = jumps_detected(ji,5);
    peakstats = get_peak_stats(Fstats,tvect,ti1,ti2,tfin); 
    npeaks = count_peaks(Fstats,ti1,ti2,Fstat_thresh_detect,Fstat_thresh_meas,0.5,1.5);   
    
    % columns are: Fstatmax, integral of F stat vs t, time of average F stat vs t,
    % time above threshold, FWHM of F stat vs t, std. dev. of F stat vs t,
    % skewness of F stat vs t, kurtosis of F stat vs t, # of peaks detected
    % in F stat vs t.
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

% plot specific jump
% [med_rel_jump, med_Fstats] = get_median_jump(~alljumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
jump_range = 1:length(med_Fstats);
ji=4;
plot_jump(nmodes,tvect(jump_range),[rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)],log10(Fstats_ts(ji,:)),log10(med_Fstats),tmeas,tjump,tjump_pre);

%% Calculate deviation from median F stat vs time

tmax = tvect(med_Fstats==max(med_Fstats));
tthresh = tvect(med_Fstats==max(med_Fstats));
% use points with time less than jump threshold crossing
tweight1 = double(tvect(1:length(med_Fstats)) < tmeas+tjump_pre)';
% use points with time less than that of the maximum F stat
tweight2 = double(tvect(1:length(med_Fstats)) < tmax)';
% columns are Euclidean distance: for the portion to the left of jump
% threshold, and for the portion to the left of maximum.
jump_stats_v_median = zeros(size(jump_stats,1),2);
for ji = 1:size(jumps_measured,1)
    
    med_Fstats_norm = log10(med_Fstats)/max(log10(med_Fstats));
    Fstats_ts_norm = log10(Fstats_ts(ji,:))/max(log10(Fstats_ts(ji,:)));
    euc_dist1 = sum(tweight1.*(Fstats_ts_norm - med_Fstats_norm).^2);
    euc_dist2 = sum(tweight2.*(Fstats_ts_norm - med_Fstats_norm).^2);
    jump_stats_v_median(ji,1) = euc_dist1;
    jump_stats_v_median(ji,2) = euc_dist2;
    
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
plotoutliers = 1;

% scatter plot of fingerprint vector
figure;
plot(jumps_measured(usejumps,2),jumps_measured(usejumps,3),'.','MarkerSize',7);
xlabel('df1'); ylabel('df2');
if plotoutliers
    hold on
    plot(jumps_measured(pickjumps1,2),jumps_measured(pickjumps1,3),'.','MarkerSize',7);
    plot(jumps_measured(pickjumps2,2),jumps_measured(pickjumps2,3),'.','MarkerSize',7);
    plot(jumps_measured(pickjumps5,2),jumps_measured(pickjumps5,3),'.','MarkerSize',7);
    plot(jumps_measured(pickjumps4,2),jumps_measured(pickjumps4,3),'.','MarkerSize',7);
end
% legend('Central cluster','Large Euclidean distance','Short FWHM','Long FWHM','Large F statistic max');

% plot fingerprint with error bars (needs to include rotation as well)
% figure; hold on
% for ji = 1:size(jumps_measured,1)
%     if ~usejumps(ji), continue; end
% %     if abs(jumps_measured(ji,4)) > .5e-6 || abs(jumps_measured(ji,5)) > .5e-6, continue; end
%     plot_ellipse(jumps_measured(ji,4),jumps_measured(ji,5),jumps_measured(ji,2),jumps_measured(ji,3),0,[0 0 1 0.2]);
% end
% xlabel('df1'); ylabel('df2');

% F stat max vs Euclidean dist
figure;
plot(jump_stats(usejumps,1),jump_stats_v_median(usejumps,1),'.');
xlabel('F statistic maximum');
ylabel('Euclidean distance');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
if plotoutliers
    hold on
    plot(jump_stats(pickjumps1,1),jump_stats_v_median(pickjumps1,1),'.');
    plot(jump_stats(pickjumps2,1),jump_stats_v_median(pickjumps2,1),'.');
    plot(jump_stats(pickjumps4,1),jump_stats_v_median(pickjumps4,1),'.');
    plot(jump_stats(pickjumps5,1),jump_stats_v_median(pickjumps5,1),'.');
end

% Euclidean distance vs FWHM
figure;
plot(jump_stats(usejumps,5),jump_stats_v_median(usejumps,1),'.','MarkerSize',7);
xlabel('FWHM (s)');
ylabel('Euclidean distance of F statistic from median');
set(gca, 'YScale', 'log');
if plotoutliers
    hold on
    plot(jump_stats(pickjumps1,5),jump_stats_v_median(pickjumps1,1),'.','MarkerSize',7);
    plot(jump_stats(pickjumps2,5),jump_stats_v_median(pickjumps2,1),'.','MarkerSize',7);
    plot(jump_stats(pickjumps5,5),jump_stats_v_median(pickjumps5,1),'.','MarkerSize',7);
    plot(jump_stats(pickjumps4,5),jump_stats_v_median(pickjumps4,1),'.','MarkerSize',7);
end

%% Manually choose outliers or jumps of interest by interacting with figures
% Use figures to manually select jumps of interest, label them as "brush1", "brush2", etc.
% then flag corresponding indices among measured jumps and confirm selection

% pick jumps on finger print scatter plot
jumps_selected1 = false(size(jumps_measured,1),1);
for ii=1:length(usejumps)
    for jj=1:size(brush1,1)
        if brush1(jj,1)==jumps_measured(ii,2) && brush1(jj,2)==jumps_measured(ii,3), jumps_selected1(ii)=1; end
    end
end

% pick jumps on Euclidean distance vs FWHM plot
jumps_selected2 = false(size(jumps_measured,1),1);
for ii=1:length(usejumps)
    for jj=1:size(brush2,1)
        if brush2(jj,1)==jump_stats(ii,5) && brush2(jj,2)==jump_stats_v_median(ii,1), jumps_selected2(ii)=1; end
    end
end

pickjumps = jumps_selected1 | jumps_selected2;

%% Manually choose outliers or jumps of interest using cutoffs of statistics

pickjumps1 = false(size(jump_stats,1),1);
pickjumps2 = false(size(jump_stats,1),1);
pickjumps3 = false(size(jump_stats,1),1);
pickjumps4 = false(size(jump_stats,1),1);
pickjumps5 = false(size(jump_stats,1),1);
pickjumps = false(size(jump_stats,1),1);

for ji = 1:length(pickjumps)
    
    Fstatmax = jump_stats(ji,1);
    t_abovethresh = jump_stats(ji,4);
    t_fwhm = jump_stats(ji,5);
    euc_dist = jump_stats_v_median(ji,1);
    euc_dist2 = jump_stats_v_median(ji,2);
    pos_jump = jumps_measured(ji,2) > 0 || jumps_measured(ji,3) > 0;
    if nmodes ==3
        pos_jump = pos_jump || jumps_measured(ji,3) > 0;
    end
    
     if euc_dist > 30
         pickjumps1(ji) = 1;
     elseif t_fwhm < .06  % .0625 slightly better?
         pickjumps2(ji) = 1;
     elseif pos_jump
         pickjumps3(ji) = 1;
     elseif Fstatmax > 7E5   % use 1E4 for 1x snr; use 7E5 for 10x snr; use 2E3 for 0.5x snr
         pickjumps4(ji) = 1;
     elseif t_fwhm > .08
         pickjumps5(ji) = 1;
     end
end

pickjumps = pickjumps1 | pickjumps2 | pickjumps3 | pickjumps4;

%% Plot median behavior of various events of interest

[med_rel_jump1, med_Fstats1] = get_median_jump(pickjumps1,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump1,log10(med_Fstats1),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Large Euclidean distance');

[med_rel_jump2, med_Fstats2] = get_median_jump(pickjumps2,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump2,med_Fstats2,med_Fstats,tmeas,tjump,tjump_pre);
title('Short FWHM');

[med_rel_jump5, med_Fstats5] = get_median_jump(pickjumps5,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump5,log10(med_Fstats5),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Long FWHM');

[med_rel_jump0, med_Fstats0] = get_median_jump(~pickjumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts);
plot_jump(nmodes,tvect(jump_range),med_rel_jump0,log10(med_Fstats0),log10(med_Fstats),tmeas,tjump,tjump_pre);
title('Non-outliers');

%% Optionally plot specific jumps of interest. TODO: needs debugging

plotjumps = find(pickjumps);%[]; % list of indexes in jumps_measured to plot

% optionally plot specific jumps of interest
for pi = 1:length(plotjumps)
    
    ji = plotjumps(pi);
    ti_jump = jumps_measured(ji,5); 
    ti = ti_jump-Nmeas;
    all_range = ti:ti+2*Nmeas+Njump-1;
    rel_jump_ts = [rel_jump_ts_1(ji,:); rel_jump_ts_2(ji,:)];
    if nmodes==3
        rel_jump_ts = [rel_jump_ts; rel_jump_ts_3(ji,:)];
    end
    plot_jump(nmodes,tvect(all_range),rel_jump_ts,Fstats(all_range),tmeas,tjump,tjump_pre);
end

%% Output chosen jumps
% have user choose output file to avoid overwriting
% writematrix(jumps_measured,uiputfile());
writematrix(jumps_measured(usejumps,:),uiputfile());

