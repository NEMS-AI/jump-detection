function peak_stats = get_peak_stats_simple(Fstats,Fstat_thresh,ti1,ti2)
%get_peak_stats summarize various statistics related to Fstat peak.
% Fstats input is the variable over entire dataset
% tvect is similarly time vector over entire dataset
% ti1: time index when jump crosses above threshold
% ti2: time index when jump crosses back below threshold

    Fstatmax = max(Fstats(ti1:ti2));
    ti_max = find(Fstats(ti1:ti2)==Fstatmax)+ti1-1;
    
    % find peak relative to right edge
    peak_i = ti2;
    while Fstats(peak_i) < (Fstatmax-Fstat_thresh)*0.5 && peak_i > ti1
        peak_i = peak_i - 1;
    end
    
    peak_stats = [Fstatmax ti_max peak_i];
end

