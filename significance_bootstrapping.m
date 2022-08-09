
noise_type = 'exp'; % 'white', 'pink', 'exp'
Nmodes = 2; % 1 or 2
p_dim = Nmodes;

m1_f0 = 4.0E7; % beginning frequency of mode 1
m2_f0 = 1.1E8;

% experimental fractional frequency noise to match to (at 100 ms pll)
m1_sigma_frac = 1.6E-7;
m2_sigma_frac = 1.0E-7;
rho = 0;  % rho = 0.19; 

tsample = .00025; % 250 microsecond sample rate
tjump = 0; % 0 or 10 ms or 50 ms for jump itself
tmeas = .1; % 100 ms to measure jump height
Nsamples = floor(tmeas/tsample);
Njump = floor(tjump/tsample);
tvect = -(Nsamples*tsample):tsample:Njump*tsample+(Nsamples-1)*tsample;
tvect = tvect';

if strcmp(noise_type,'white')
    % to find adev at tsample seconds given sigma_frac adev at 0.1 seconds:
    % a = .3162278 * sigma_frac * tsample^(-.5)
    adev1 = .3162278 * m1_sigma_frac * tsample^(-.5);
    adev2 = .3162278 * m2_sigma_frac * tsample^(-.5);
else
    adev1 = m1_sigma_frac;
    adev2 = m2_sigma_frac;
    if strcmp(noise_type,'pink')
        n1 = importdata('n1_pink.txt');
        n2 = importdata('n2_pink.txt');
    else
        n1 = importdata('n1_detrend.txt');
        n2 = importdata('n2_detrend.txt');
    end
    n1_frac = n1(:,2)/mean(n1(:,2))-1;
    n2_frac = n2(:,2)/mean(n2(:,2))-1;
end

% frequency noise in Hz
sigma1 = adev1*m1_f0;
sigma2 = adev2*m2_f0;
% noise covariance matrix (over entire dataset), used only to generate white noise
SIGMA = [sigma1^2, rho*sigma1*sigma2;
         rho*sigma1*sigma2, sigma2^2];

% optionally plot fractional frequency
% figure;
% plot(n1(:,1),n1_frac);
% hold on
% plot(n2(:,1),n2_frac);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');

% fraction under x sigma (analog for normal distribution). And area to look for 2 tail test.
frac_1sig = 0.682689492137086; frac_1sig_2tail = 1-(1-frac_1sig)/2;
frac_2sig = 0.954499736103642; frac_2sig_2tail = 1-(1-frac_2sig)/2;
frac_3sig = 0.997300203936740; frac_3sig_2tail = 1-(1-frac_3sig)/2;

Niter = 10000;
stats = zeros(Niter,1);

sigma_samples = [];

for i=1:Niter
         
    if strcmp(noise_type,'white')
        nvect = mvnrnd([0,0],SIGMA,length(tvect));
        f1 = m1_f0 + nvect(:,1);
        f2 = m2_f0 + nvect(:,2);
    else
        % randomly sample (with replacement) from noise data
        ti_start = floor(rand()*(size(n1,1)-2*Nsamples-Njump-100));
        f1 = m1_f0*(1+n1_frac(ti_start:ti_start+length(tvect)-1));
        f2 = m2_f0*(1+n2_frac(ti_start:ti_start+length(tvect)-1));
    end
    
    if Nmodes == 1
        xi = f1(1:Nsamples); yi = f1(Njump+Nsamples+1:end);
        sigma_pool_sample = sqrt((var(xi)+var(yi))/2);
        sigma_samples = [sigma_samples; sigma_pool_sample];
        t_m1 = (mean(xi)-mean(yi)) / (sigma_pool_sample * sqrt(2/Nsamples));
        stats(i) = t_m1;
    else
        xi = [f1(1:Nsamples)'; f2(1:Nsamples)'];
        yi = [f1(Njump+Nsamples+1:end)'; f2(Njump+Nsamples+1:end)'];
        xbar = mean(xi,2);
        ybar = mean(yi,2);
        sigma_x = cov(xi');
        sigma_y = cov(yi');
        sigma_pool = sigma_x/2 + sigma_y/2;
        sigma_samples = [sigma_samples; sigma_pool(1,1) sigma_pool(1,2) sigma_pool(2,2)];
        t2 = Nsamples/2*(xbar-ybar)'/sigma_pool*(xbar-ybar);
        Fstat = (2*Nsamples-p_dim-1)/(p_dim*(2*Nsamples-2))*t2;
        stats(i) = Fstat;
    end
end

stats_sort = sort(stats);
stats_frac = linspace(0,1,Niter);
figure;
plot(stats_sort,stats_frac);
if Nmodes == 1
    stat_1sig = stats_sort(find(stats_frac > frac_1sig_2tail,1));
    stat_2sig = stats_sort(find(stats_frac > frac_2sig_2tail,1));
    stat_3sig = stats_sort(find(stats_frac > frac_3sig_2tail,1));
else
    stat_1sig = stats_sort(find(stats_frac > frac_1sig,1));
    stat_2sig = stats_sort(find(stats_frac > frac_2sig,1));
    stat_3sig = stats_sort(find(stats_frac > frac_3sig,1));
end

% calculate frequency shift corresponding to 1 sigma significance
if Nmodes == 1
    tinv_1sig = stat_1sig;
    fmeandiff_1sig = tinv_1sig * mean(sigma_samples) * sqrt(2/Nsamples);
    fmeandiff_1sig / m1_f0
else
    t2inv_1sig = stat_1sig/(2*Nsamples-p_dim-1)*(p_dim*(2*Nsamples-2));
%     t2inv_1sig = stat_1sig;

    % in 2 dimensions we get an ellipse.
    S11 = mean(sigma_samples(:,1));
    S12 = mean(sigma_samples(:,2));
    S22 = mean(sigma_samples(:,3));
    A = S22;
    B = -2*S12;
    C = S11;
    F = -t2inv_1sig  / Nsamples * (S11*S22-S12^2);

    % semi-major axis
    a = -sqrt(2*((B^2 - 4*A*C)*F)*((A+C) + sqrt((A-C)^2 + B^2)))/(B^2-4*A*C); % semi-major
    % semi-minor axis
    b = -sqrt(2*((B^2 - 4*A*C)*F)*((A+C) - sqrt((A-C)^2 + B^2)))/(B^2-4*A*C); % semi-minor
    % the angle from the positive horizontal axis to the ellipse's major axis
    if B ~= 0
        theta = atan(1/B*(C-A-sqrt((A-C)^2 + B^2)));
    elseif A < C
        theta = 0;
    else
        theta = pi/2;
    end
    
    fmeandiff_1sig = zeros(Niter,2);
    for i = 1:Niter
        theta_i = (i-1)/Niter*2*pi;
        fmeandiff_1sig(i,:) = [a*cos(theta_i)*cos(theta)-b*sin(theta_i)*sin(theta) ...
                               a*cos(theta_i)*sin(theta)-b*sin(theta_i)*cos(theta)];
    end

    figure;
    plot(fmeandiff_1sig(:,1)/m1_f0,fmeandiff_1sig(:,2)/m2_f0);
    sqrt(pi*a*b/(m1_f0*m2_f0))
end

signif_1sig = zeros(Niter,1);

for i=1:Niter
    
    if strcmp(noise_type,'white')
        nvect = mvnrnd([0,0],SIGMA,length(tvect));
        f1 = m1_f0 + nvect(:,1);
        f2 = m2_f0 + nvect(:,2);
    else
        % randomly sample (with replacement) from noise data
        ti_start = floor(rand()*(size(n1,1)-2*Nsamples-Njump-100));
        f1 = m1_f0*(1+n1_frac(ti_start:ti_start+length(tvect)-1));
        f2 = m2_f0*(1+n2_frac(ti_start:ti_start+length(tvect)-1));
    end
    
    % jump size
    if Nmodes == 1, f1_shift = fmeandiff_1sig;
    else, f1_shift = fmeandiff_1sig(i,1); end
    
    f1_syn = [m1_f0*ones(Nsamples,1); ...
             (m1_f0-f1_shift)*ones(Njump,1); ...
             (m1_f0-f1_shift)*ones(Nsamples,1)]; 
    f1 = f1 + f1_syn;
    
    if Nmodes == 2
        f2_shift = fmeandiff_1sig(i,2);
        f2_syn = [m2_f0*ones(Nsamples,1); ...
                 (m2_f0-f2_shift)*ones(Njump,1); ...
                 (m2_f0-f2_shift)*ones(Nsamples,1)];
        f2 = f2 + f2_syn;
    end
    
    if Nmodes == 1
        xi = f1(1:Nsamples); yi = f1(Njump+Nsamples+1:end);
        sigma_pool_sample = sqrt((var(xi)+var(yi))/2);
        t_m1 = (mean(xi)-mean(yi)) / (sigma_pool_sample * sqrt(2/Nsamples));
        signif_1sig(i) = t_m1 > stat_1sig;
    else
        xi = [f1(1:Nsamples)'; f2(1:Nsamples)'];
        yi = [f1(Njump+Nsamples+1:end)'; f2(Njump+Nsamples+1:end)'];
        xbar = mean(xi,2);
        ybar = mean(yi,2);
        sigma_x = cov(xi');
        sigma_y = cov(yi');
        sigma_pool = sigma_x/2 + sigma_y/2;
        t2 = Nsamples/2*(xbar-ybar)'/sigma_pool*(xbar-ybar);
        Fstat = (2*Nsamples-p_dim-1)/(p_dim*(2*Nsamples-2))*t2;
        signif_1sig(i) = Fstat > stat_1sig;
    end
end

sum(signif_1sig)/Niter
