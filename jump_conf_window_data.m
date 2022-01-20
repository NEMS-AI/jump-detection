clear all
% fraction under x sigma. And area to look for 2 tail test.
% frac_1sig = 0.682689492137086; frac_1sig_2tail = 1-(1-frac_1sig)/2;
% frac_2sig = 0.954499736103642; frac_2sig_2tail = 1-(1-frac_2sig)/2;
% frac_3sig = 0.997300203936740; frac_3sig_2tail = 1-(1-frac_3sig)/2;
% frac_4sig = 0.999936657516334; frac_4sig_2tail = 1-(1-frac_4sig)/2;
% frac_5sig = 0.999999426696856; frac_5sig_2tail = 1-(1-frac_5sig)/2;

% % read noise
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

% read data
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
data = readtable('synthetic/2/5_inst_exp_data.csv');
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

plot(tvect,f1/f1(1),'k'); hold on
plot(tvect,f2/f2(1),'b');

fdiff = [];

tsample = .00025; % 250 microsecond sample rate
tjump = .0; % 0 or 10 ms or 50 ms for jump itself
tmeas = .050; % 100 ms to measure jump height
Nsamples = floor(tmeas/tsample);
Njump = floor(tjump/tsample);

tfin = ttot-Nsamples*2-Njump;

Fstats = zeros(tfin,1);
t_Fstats = tvect(Nsamples:tfin+Nsamples-1);
pvals = [];

fprintf('jump_conf:         ');
ti = 1;

jump_events = [];

while ti < tfin
   
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
    if Fstat < 1200
        ti = ti+1;
        Fstats(ti) = Fstat;
    else
        % find max Fstat within next Nsamples
        Fstatmax = 0;
        ti_max = 0;
        rel_jump = 0;
        for tii = 0:Nsamples-1
            xii_range = xi_range+tii;
            yii_range = yi_range+tii;

            fvect_xii = fvect(:,xii_range);
            fvect_yii = fvect(:,yii_range);

            xbar = mean(fvect_xii,2);
            ybar = mean(fvect_yii,2);
            sigma_x = cov(fvect_xii');
            sigma_y = cov(fvect_yii');
            sigma_pool = sigma_x/2 + sigma_y/2;
            t2 = Nsamples/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
            p_dim = 2;
            Fstat_ii = (2*Nsamples-p_dim-1)/(p_dim*(2*Nsamples-2))*t2;       
            % sometimes the below causes len(Fstats) > tfin
            Fstats(ti+tii) = Fstat_ii;
            if Fstat_ii > Fstatmax
                Fstatmax = Fstat_ii;
                ti_max = ti+tii+Nsamples-1; %-1: event happens 1 step before detect
                rel_jump = ybar./xbar-1;
            end
        end
        jump_events = [jump_events; ti_max tvect(ti_max) rel_jump' Fstatmax];
        1;
        ti = ti + Nsamples + 1;
    end
    % not sure if below is correct.
%     pval = 1 - fcdf(Fstat, p_dim, 2*Nsamples-p_dim-1);
    
end
%%
% redo calculation for jumps that are tmeas2 apart.
tmeas2 = .5; % 500 ms to measure jump height
Nsamples2 = floor(tmeas2/tsample);
jump_events2 = [];

for ji = 3:size(jump_events,1)-3
        
    t_prev = tvect(jump_events(ji-1,1));
    t_curr = tvect(jump_events(ji,1));
    t_next = tvect(jump_events(ji+1,1));
    
    if t_curr-t_prev < tmeas2 || t_next-t_curr < tmeas2
        continue
    end
    
    ti_center = jump_events(ji,1);
    ti = ti_center-Nsamples2;
    xi_range = ti:ti+Nsamples2-1;
    yi_range = ti+Nsamples2+Njump:ti+2*Nsamples2+Njump-1;

    fvect_xi = fvect(:,xi_range);
    fvect_yi = fvect(:,yi_range);

    xbar = mean(fvect_xi,2);
    ybar = mean(fvect_yi,2);
    
    rel_jump = ybar./xbar-1;
    jump_events2 = [jump_events2; tvect(ti_center) rel_jump'];
end
% yyaxis right
% plot(t_Fstats,Fstats(1:end));
