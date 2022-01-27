clear all

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
% % scale noise
% f1_start = f1(1); f2_start = f2(2);
% f1 = (f1-f1_start)*noise_factor + f1_start;
% f2 = (f2-f2_start)*noise_factor + f2_start;



% % read data
% fid = fopen('experiment/[CL m1] 13_Z=4665 X=2388 Y=3090 #20 GroEl.txt');
% nems1 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% fid = fopen('experiment/[CL m2] 13_Z=4665 X=2388 Y=3090 #20 GroEl.txt');
% nems2 = textscan(fid, '%f %f %f %f %f','HeaderLines', 82 );
% fclose(fid);
% tvect=nems1{1}; f1=nems1{4};
% tvect=tvect-tvect(1);
% %t2=nems2{1};
% f2=nems2{4};

% read synthetic data
data = readtable('synthetic/2/7_decay_exp_data.csv');
data = table2array(data);
tvect = data(:,1); tvect=tvect-tvect(1);
f1 = data(:,2);
f2 = data(:,3);

% events = readtable('synthetic/2/inst_exp_events.csv');
% events = table2array(events);
% events(:,2) = events(:,2)/f1(1);
% events(:,3) = events(:,3)/f2(1);

% preprocess to remove drift? invert the below:
% f1_drift = (.99997-1)*mode1_fstart/700*tvect;
% f2_drift = (.99998-1)*mode2_fstart/500*tvect;

ttot = length(tvect);
tvect=tvect(1:ttot);
f1=f1(1:ttot);
f2=f2(1:ttot);
fvect=[f1 f2]';
clear nems1 nems2

% plot(tvect,f1/f1(1),'k'); hold on
% plot(tvect,f2/f2(1),'b');

fdiff = [];

tsample = .00025; % 250 microsecond sample rate
tjump = .04; % 0 or 10 ms or 50 ms for jump itself
tmeas = .10; % 100 ms to measure jump height
Nsamples = floor(tmeas/tsample);
Njump = floor(tjump/tsample);

tfin = ttot-Nsamples*3-Njump;
% tfin = floor(tfin/10);

Fstats = zeros(tfin,1);
t_Fstats = tvect(Nsamples:tfin+Nsamples-1);
pvals = [];
Fstat_thresh = 2e3;

fprintf('Fstat:         ');
ti = 1;
for ti=1:tfin
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 

    xi_range = ti:ti+Nsamples-1;
    yi_range = ti+Nsamples+Njump:ti+2*Nsamples+Njump-1;

    fvect_xi = fvect(:,xi_range);
    fvect_yi = fvect(:,yi_range);

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    sigma_x = cov(fvect_xi');
    sigma_y = cov(fvect_yi');
    sigma_pool = sigma_x/2 + sigma_y/2;
    t2 = Nsamples/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
    p_dim = 2;
    Fstat = (2*Nsamples-p_dim-1)/(p_dim*(2*Nsamples-2))*t2;
    Fstats(ti) = Fstat;
    
end
%%
fprintf('jump_events:         ');
jump_events = [];
ti=1;
while ti < tfin
    
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/ttot*100); 
    if Fstats(ti) < Fstat_thresh
        ti = ti+1;
    else
        tii = 0;
        while tii < tfin-ti && Fstats(ti+tii) > Fstat_thresh
            tii = tii+1;
        end
        peak_width = tii;
        Fi_max = floor(ti+tii/2);
        Fstatmax = Fstats(Fi_max);
        
        % narrow to FWHM location
        peak_left_i = ti;
        while Fstats(peak_left_i) < Fstatmax/2 && peak_left_i < tfin
            peak_left_i = peak_left_i + 1;
        end
        peak_right_i = ti+tii;
        while Fstats(peak_right_i) < Fstatmax/2 && peak_right_i > 1
            peak_right_i = peak_right_i - 1;
        end
        
        peak_width = peak_right_i - peak_left_i;
%         if tvect(peak_width) > tjump*.75
%             ti = ti + tii + 1;
%             continue;
%         end
%         Fi_max = floor((peak_left_i+peak_right_i)/2);
%         Fi_max = floor(peak_right_i);
        ti_max = ti+Nsamples-1+tii;
        jump_events = [jump_events; ti_max tvect(ti_max) tvect(peak_width) Fstatmax];
        ti = ti + tii + 1;
    end
end

%%
% do final calculation for jumps that are tmeas2 apart.
tmeas2 = .2; % 500 ms to measure jump height
Nsamples2 = floor(tmeas2/tsample);
jump_events2 = [];
rel_jump_ts_1 = [];
rel_jump_ts_2 = [];

for ji = 2:size(jump_events,1)-2
        
    t_prev = tvect(jump_events(ji-1,1));
    t_curr = tvect(jump_events(ji,1));
    t_next = tvect(jump_events(ji+1,1));
    
    if t_curr-t_prev < tmeas2 || t_next-t_curr < tmeas2, continue; end
    
    ti_center = jump_events(ji,1);
    ti = ti_center-Nsamples2;
    xi_range = ti:ti+Nsamples2-1;
    yi_range = ti+Nsamples2+Njump:ti+2*Nsamples2+Njump-1;
    all_range = ti:ti+2*Nsamples2+Njump-1;

    fvect_xi = fvect(:,xi_range);
    fvect_yi = fvect(:,yi_range);

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    
    rel_jump = ybar./xbar-1;
    jump_events2 = [jump_events2; tvect(ti_center) rel_jump'];
    
    rel_jump_ts = (fvect(:,all_range)-ybar)./(xbar-ybar);
    rel_jump_ts_1 = [rel_jump_ts_1; rel_jump_ts(1,:)];
    rel_jump_ts_2 = [rel_jump_ts_2; rel_jump_ts(2,:)];
    
%     figure;
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas2,rel_jump_ts(1,:),'k'); hold on
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas2,rel_jump_ts(2,:),'b');
%     yyaxis right
%     plot(tvect(1:size(rel_jump_ts,2))-tmeas2,Fstats(all_range-Nsamples));
%     1;
end

med_rel_jump = [median(rel_jump_ts_1); median(rel_jump_ts_2);];

%%
% avg_rel_jump = sum_rel_jump / size(jump_events,1);
figure;
xtofit=tvect(1:size(med_rel_jump,2))-tmeas2;
ytofit1=med_rel_jump(1,:);
ytofit2=med_rel_jump(2,:);
xtofit_aj = xtofit(Nsamples+1:end);
ytofit1_aj = ytofit1(Nsamples+1:end);
ytofit2_aj = ytofit2(Nsamples+1:end);
plot(xtofit,ytofit1,'k'); hold on
plot(xtofit,ytofit2,'b');
% experimental fit
% plot(xtofit_aj,.946*exp(-xtofit_aj/.01544),'r');
% synthetic fit, 1/10 and 1/100
% plot(xtofit_aj,.9556*exp(-101.9*xtofit_aj),'r');
% synthetic fit
% plot(xtofit_aj,.949*exp(-62.51*xtofit_aj),'r');
% plot(xtofit_aj,.7312*exp(-104.4*xtofit_aj),'r');
1;

