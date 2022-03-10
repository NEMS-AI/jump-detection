f = figure;
axis image;


%% All plots utilize the following
% scatter plot of fingerprint vector
subplot(2,2,1);
plot(jumps_measured(usejumps,2),jumps_measured(usejumps,3),'.');
xlabel('df1'); ylabel('df2');
if plotoutliers
    hold on
    plot(jumps_measured(pickjumps,2),jumps_measured(pickjumps,3),'.');
end

% F stat max vs time above threshold
subplot(2,2,2);
plot(jump_stats(usejumps,1),jump_stats_v_median(usejumps,1),'.');
xlabel('F statistic maximum');
ylabel('Euclidean distance');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
if plotoutliers
    hold on
    plot(jump_stats(pickjumps,1),jump_stats_v_median(pickjumps,1),'.');
end

% FWHM vs time above threshold
subplot(2,2,3);
plot(jump_stats(usejumps,4),jump_stats(usejumps,5),'.');
xlabel('Time above threshold (s)');
ylabel('FWHM (s)');
if plotoutliers
    hold on
    plot(jump_stats(pickjumps,4),jump_stats(pickjumps,5),'.');
end


