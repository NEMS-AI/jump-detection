
nmodes = 2;
data_type = 'noise';
tsample = .00025;        % seconds per sample
tmeas_detect = .10;      % time window for detecting jump (before and after samples)
tjump = .06;             % time window for (full) jump itself
tjump_pre = .03;         % time before F threshold detection for jump
tjump_post = tjump - tjump_pre; % time after detection for jump itself 

% number of samples for detection
Ndetect = floor(tmeas_detect/tsample); 
Njump = floor(tjump/tsample);
Npre = floor(tjump_pre/tsample);

% % read noise
% noise_factor = .01;
% fid = fopen('experiment/[CL m1] 12_Z=4665 X=2388 Y=3090 #20 no ions.txt');
% nems1 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% fid = fopen('experiment/[CL m2] 12_Z=4665 X=2388 Y=3090 #20 no ions.txt');
% nems2 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% tvect=nems1{1}; f1=nems1{4};
% tvect=tvect-tvect(1);
% %t2=nems2{1};
% f2=nems2{4};
% 
% % scale noise for simulation purpose
% f1_start = f1(1); f2_start = f2(2);
% f1 = (f1-f1_start)*noise_factor + f1_start;
% f2 = (f2-f2_start)*noise_factor + f2_start;

[tvect, fvect] = load_data(filenames, data_type, nmodes);



% % read data
% outfile = 'experiment/GroEL events.csv';
% fid = fopen('experiment/[CL m1] 13_Z=4665 X=2388 Y=3090 #20 GroEl.txt');
% nems1 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% fid = fopen('experiment/[CL m2] 13_Z=4665 X=2388 Y=3090 #20 GroEl.txt');
% nems2 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% tvect=nems1{1}; f1=nems1{4};
% tvect=tvect-tvect(1);
% % tvect2=nems2{1}; 
% % tvect2=tvect2-tvect2(1);
% f2=nems2{4};
% tvectuse = tvect > 10 & tvect < 830;  % full dataset
% % tvectuse = tvect > 10 & tvect < 420;  % first half
% % tvectuse = tvect > 420 & tvect < 830;  % second half
% % tvectuse = tvect > 10 & tvect < 92;  % first tenth
% % tvectuse = tvect > 748 & tvect < 830;  % last tenth
% tvect = tvect(tvectuse);
% tvect = tvect-tvect(1);
tstart = 10; % used for time lag calculation
% f1=f1(tvectuse);
% f2=f2(tvectuse);
tlag_mode2_per_second = .013/830;

% % read synthetic data
% infile = 'synthetic/2/7_decay_exp_data.csv';
% outfile = 'synthetic/2/7_decay_exp_detected.csv';
% data = readtable(infile);
% data = table2array(data);
% tvect = data(:,1); tvect=tvect-tvect(1);
% f1 = data(:,2);
% f2 = data(:,3);

% preprocess to remove drift? invert the below:
% f1_drift = (.99997-1)*mode1_fstart/700*tvect;
% f2_drift = (.99998-1)*mode2_fstart/500*tvect;

ttot = length(tvect);

% last time point
tfin = ttot-Ndetect*3-Njump;
% tfin = floor(tfin/50);   % optionally use subset of data to peek at results

Fstats = zeros(tfin,1);
t_Fstats = tvect(Ndetect:tfin+Ndetect-1);
Fstat_thresh_detect = 1200;
Fstat_thresh_meas = 2400;

% plot(tvect,f1/f1(1),'k'); hold on
% plot(tvect,f2/f2(1),'b');

%%

fprintf('Fstat:         ');
for ti=100:tfin-100 % allow for time lag in either direction
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 

    xi_range = ti:ti+Ndetect-1;
    yi_range = ti+Ndetect+Njump:ti+2*Ndetect+Njump-1;
    
    tlag_m2 = (tvect(ti)+tstart)*tlag_mode2_per_second;
    Nlag_m2 = floor(tlag_m2/tsample);

    fvect_xi = [fvect(1,xi_range); fvect(2,xi_range+Nlag_m2)];
    fvect_yi = [fvect(1,yi_range); fvect(2,yi_range+Nlag_m2)];

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    sigma_x = cov(fvect_xi');
    sigma_y = cov(fvect_yi');
    sigma_pool = sigma_x/2 + sigma_y/2;
    t2 = Ndetect/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
    p_dim = 2;
    Fstat = (2*Ndetect-p_dim-1)/(p_dim*(2*Ndetect-2))*t2;
    Fstats(ti+Ndetect) = Fstat;
    
end
%%
fprintf('jump_events:         ');
jumps_detected = [];
ti=100; % allow for time lag in either direction
while ti < tfin-100
    
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
tmeas = .2; % 200 ms to measure jump height
Nmeas = floor(tmeas/tsample);
jumps_measured = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];
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

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    
    rel_jump = ybar./xbar-1;
    jumps_measured = [jumps_measured; tvect(ti_jump) peak_width_fwhm Fstatmax rel_jump'];
    
    rel_jump_ts = ([fvect(1,all_range); fvect(2,all_range+Nlag_m2);]-ybar)./(xbar-ybar);
    rel_jump_ts_1 = [rel_jump_ts_1; rel_jump_ts(1,:)];
    rel_jump_ts_2 = [rel_jump_ts_2; rel_jump_ts(2,:)];
    Fstats_ts = [Fstats_ts; Fstats(all_range)'];
    
%     figure;
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,rel_jump_ts(1,:),'k'); hold on
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,rel_jump_ts(2,:),'b');
%     yyaxis right
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas-tjump_pre,Fstats(all_range));%-Nsamples));
%     tvect_ts = tvect(all_range);
%     1;
end

writematrix(jumps_measured,outfile);

%%
% plot jump signature with optional exponential fit (obtained separately with toolbox)
med_rel_jump = [median(rel_jump_ts_1); median(rel_jump_ts_2);];
med_Fstats = median(Fstats_ts);

figure;
tvect_med_jump=tvect(1:size(med_rel_jump,2))-tmeas;
fvect_med_jump_1=med_rel_jump(1,:);
fvect_med_jump_2=med_rel_jump(2,:);
plot(tvect_med_jump,fvect_med_jump_1,'k'); hold on
plot(tvect_med_jump,fvect_med_jump_2,'b');
xlabel('Time (s)');
ylabel('Normalized jump');
legend('Mode 1','Mode 2');
% yyaxis right
% plot(tvect_med_jump,med_Fstats);

% portion after the jump to fit exponential to
t0_exp = .0125;
xtofit_aj = tvect_med_jump(tvect_med_jump > t0_exp)- t0_exp;
ytofit1_aj = fvect_med_jump_1(tvect_med_jump > t0_exp);
ytofit2_aj = fvect_med_jump_2(tvect_med_jump > t0_exp);
% synthetic fit, 1/10 and 1/100
% plot(xtofit_aj+t0_exp,1.008*exp(-99.83*xtofit_aj),'r');
% synthetic fit, 1
% plot(xtofit_aj+t0_exp,.9037*exp(-90.33*xtofit_aj),'r');
% experimental fit
plot(xtofit_aj+t0_exp,.8874*exp(-138.1*xtofit_aj),'r');
1;

