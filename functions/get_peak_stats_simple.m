function peak_stats = get_peak_stats_simple(Fstats,ti1,ti2,tfin)
%get_peak_stats summarize various statistics related to Fstat peak.
% Fstats input is the variable over entire dataset
% tvect is similarly time vector over entire dataset
% ti1: time index when jump crosses threshold
% ti2: time index when jump crosses back below threshold
% (we assume Fstat is above threshold between ti1 and ti2)

    Fstatmax = max(Fstats(ti1:ti2));
    ti_max = find(Fstats(ti1:ti2)==Fstatmax)+ti1-1;
    
    %  find left edge 
    peak_left_i = ti_max;
    while Fstats(peak_left_i) > Fstatmax*0.1 && peak_left_i > 1
        peak_left_i = peak_left_i - 1;
    end

    % find peak relative to left edge
    peak_i1 = peak_left_i;
    while Fstats(peak_i1) < Fstatmax*0.9 && peak_i1 < tfin
        peak_i1 = peak_i1 + 1;
    end
    
    %  find right edge 
    peak_right_i = ti_max;
    while Fstats(peak_right_i) > Fstatmax*0.1 && peak_right_i < tfin
        peak_right_i = peak_right_i + 1;
    end

    % find peak relative to right edge
    peak_i2 = peak_right_i;
    while Fstats(peak_i2) < Fstatmax*0.9 && peak_i2 > 1
        peak_i2 = peak_i2 - 1;
    end
    
    peak_stats = [Fstatmax ti_max peak_i1 peak_i2];
end

