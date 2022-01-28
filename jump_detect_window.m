
nmodes = 3;
data_type = 'data';     % choose from 'noise', 'data', 'synthetic'
if nmodes == 2
    tsample = .00025;        % seconds per sample
    tmeas_detect = .10;      % time window for detecting jump (before and after samples)
    tjump = .06;             % time window for (full) jump itself
    tjump_pre = .03;         % time before F threshold detection for jump
    tmeas = .2;              % 200 ms for final measurement of jump height
elseif nmodes == 3
    tsample = 0.02;
    tmeas_detect = 0.5;     % 1-10s optimal detection window
    tjump = 0.1;            % not sure exact PLL settings
    tjump_pre = 0.05;
    tmeas = 2; 
end
tjump_post = tjump - tjump_pre; % time after detection for jump itself 

% number of samples for detection
Ndetect = floor(tmeas_detect/tsample); 
Njump = floor(tjump/tsample);
Npre = floor(tjump_pre/tsample);

[tvect, fvect] = load_data(filenames, data_type, nmodes);

% optionally select portion of data
if strcmp(data_type, 'data')
    if nmodes == 2
        % tvectuse = tvect > 10 & tvect < 830;  % full dataset
        % tvectuse = tvect > 10 & tvect < 420;  % first half
        % tvectuse = tvect > 420 & tvect < 830;  % second half
        tvectuse = tvect > 10 & tvect < 92;  % first tenth
        % tvectuse = tvect > 748 & tvect < 830;  % last tenth
        tvect = tvect(tvectuse);
        tvect = tvect-tvect(1);
        fvect = fvect(:,tvectuse);
        tstart = 10; % used for time lag calculation
        tlag_mode2_per_second = .013/830;
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
end

% preprocess to remove drift? invert the below:
% f1_drift = (.99997-1)*mode1_fstart/700*tvect;
% f2_drift = (.99998-1)*mode2_fstart/500*tvect;

ttot = length(tvect);

% last time point
tfin = ttot-Ndetect*3-Njump;
% tfin = floor(tfin/50);   % optionally use subset of data to peek at results

Fstats = zeros(tfin,1);
% Fstat_thresh_detect = 1200;
% Fstat_thresh_meas = 2400;
Fstat_thresh_detect = 300;
Fstat_thresh_meas = 600;

% optionally plot relative frequencies
% figure;
% plot(tvect,fvect(1,:)/fvect(1,1),'k'); hold on
% plot(tvect,fvect(2,:)/fvect(2,1),'b');
% plot(tvect,fvect(3,:)/fvect(3,1),'r');
% xlabel('Time (s)')
% ylabel('Relative frequency');

%%

fprintf('Fstat:         ');
for ti=1+tlag_buffer:tfin-tlag_buffer % allow for time lag in either direction
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 

    xi_range = ti:ti+Ndetect-1;
    yi_range = ti+Ndetect+Njump:ti+2*Ndetect+Njump-1;
    
    tlag_m2 = (tvect(ti)+tstart)*tlag_mode2_per_second;
    Nlag_m2 = floor(tlag_m2/tsample);

    fvect_xi = [fvect(1,xi_range); fvect(2,xi_range+Nlag_m2)];
    fvect_yi = [fvect(1,yi_range); fvect(2,yi_range+Nlag_m2)];
    
    if nmodes == 3
        tlag_m3 = (tvect(ti)+tstart)*tlag_mode3_per_second;
        Nlag_m3 = floor(tlag_m3/tsample);
        fvect_xi = [fvect_xi; fvect(3,xi_range+Nlag_m3)];
        fvect_yi = [fvect_yi; fvect(3,yi_range+Nlag_m3)];
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

% optionally plot Fstats
figure;
plot(tvect(1:length(Fstats)),Fstats);

%%
fprintf('jump_events:         ');
jumps_detected = [];
ti=1+tlag_buffer; % allow for time lag in either direction
while ti < tfin-tlag_buffer
    
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 
    if Fstats(ti) < Fstat_thresh_detect
        ti = ti+1;
    else
        tii = 0;
        while tii < tfin-ti && Fstats(ti+tii) > Fstat_thresh_detect
            tii = tii+1;
        end
        Fstatmax = max(Fstats(ti:ti+tii));
        
        % find approx FWHM
        peak_left_i = ti;
        while Fstats(peak_left_i) < Fstatmax/2 && peak_left_i < tfin
            peak_left_i = peak_left_i + 1;
        end
        peak_right_i = ti+tii;
        while Fstats(peak_right_i) < Fstatmax/2 && peak_right_i > 1
            peak_right_i = peak_right_i - 1;
        end
        
        peak_width_fwhm_range = peak_right_i - peak_left_i;
        % -1 is convention: jump happens 1 sample before it can be measured
        ti_jump = ti+tii-1;
        jumps_detected = [jumps_detected; ti_jump tvect(ti_jump) tvect(peak_width_fwhm_range) Fstatmax];
        ti = ti + tii + 1;
    end
end

%%
% do final calculation for jumps that are tmeas apart.
Nmeas = floor(tmeas/tsample);
jumps_measured = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];
rel_jump_ts_3 = [];
Fstats_ts = [];

for ji = 2:size(jumps_detected,1)-2
        
    t_prev = tvect(jumps_detected(ji-1,1));
    t_curr = tvect(jumps_detected(ji,1));
    t_next = tvect(jumps_detected(ji+1,1));
    
    Fstatmax = jumps_detected(ji,4);
    if Fstatmax < Fstat_thresh_meas, continue; end
    
    if t_curr - tjump_pre < t_prev + tmeas + tjump_post || ...
       t_next - tjump_pre < t_curr + tmeas + tjump_post, continue; end
    
    % can optionally exclude jumps with too wide peak width
    peak_width_fwhm = jumps_detected(ji,3);
    if peak_width_fwhm > tjump*1.5, continue; end
    
    ti_jump = jumps_detected(ji,1);
    ti = ti_jump-Npre-Nmeas;
    xi_range = ti:ti+Nmeas-1;
    yi_range = ti+Nmeas+Njump:ti+2*Nmeas+Njump-1;
    all_range = ti:ti+2*Nmeas+Njump-1;
    
    tlag_m2 = (tvect(ti)+tstart)*tlag_mode2_per_second;
    Nlag_m2 = floor(tlag_m2/tsample);

    fvect_xi = [fvect(1,xi_range); fvect(2,xi_range+Nlag_m2)];
    fvect_yi = [fvect(1,yi_range); fvect(2,yi_range+Nlag_m2)];
    
    if nmodes == 3
        tlag_m3 = (tvect(ti)+tstart)*tlag_mode3_per_second;
        Nlag_m3 = floor(tlag_m3/tsample);
        fvect_xi = [fvect_xi; fvect(3,xi_range+Nlag_m3)];
        fvect_yi = [fvect_yi; fvect(3,yi_range+Nlag_m3)];
    end

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    
    rel_jump = ybar./xbar-1;
    jumps_measured = [jumps_measured; tvect(ti_jump) peak_width_fwhm Fstatmax rel_jump'];
    
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
    
%     figure;
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,rel_jump_ts(1,:),'k'); hold on
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,rel_jump_ts(2,:),'b');
%     yyaxis right
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,Fstats(all_range));%-Nsamples));
%     tvect_ts = tvect(all_range);
%     1;
end

% have user choose output file to avoid overwriting
writematrix(jumps_measured,uiputfile());

%%
% plot jump signature with optional exponential fit (obtained separately with toolbox)
med_rel_jump = [median(rel_jump_ts_1); median(rel_jump_ts_2)];
if nmodes == 3
    med_rel_jump = [med_rel_jump; median(rel_jump_ts_3)];
end
med_Fstats = median(Fstats_ts);

figure;
tvect_med_jump=tvect(1:size(med_rel_jump,2))-tmeas;
fvect_med_jump_1=med_rel_jump(1,:);
fvect_med_jump_2=med_rel_jump(2,:);
plot(tvect_med_jump,fvect_med_jump_1,'k'); hold on
plot(tvect_med_jump,fvect_med_jump_2,'b');
if nmodes == 3
    fvect_med_jump_3=med_rel_jump(3,:);
    plot(tvect_med_jump,fvect_med_jump_3,'r');
end
xlabel('Time (s)');
ylabel('Normalized jump');
legend('Mode 1','Mode 2','Mode 3');
% yyaxis right
% plot(tvect_med_jump,med_Fstats);

% % portion after the jump to fit exponential to
% t0_exp = .0125;
% xtofit_aj = tvect_med_jump(tvect_med_jump > t0_exp)- t0_exp;
% ytofit1_aj = fvect_med_jump_1(tvect_med_jump > t0_exp);
% ytofit2_aj = fvect_med_jump_2(tvect_med_jump > t0_exp);
% % synthetic fit, 1/10 and 1/100
% % plot(xtofit_aj+t0_exp,1.008*exp(-99.83*xtofit_aj),'r');
% % synthetic fit, 1
% % plot(xtofit_aj+t0_exp,.9037*exp(-90.33*xtofit_aj),'r');
% % experimental fit
% plot(xtofit_aj+t0_exp,.8874*exp(-138.1*xtofit_aj),'r');
1;

