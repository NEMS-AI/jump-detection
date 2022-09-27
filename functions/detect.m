function jumps_detected = detect(params, Fstats)

jumps_detected = [];

% Iterate iterate through the sequence of fstats
ti=1; 
fprintf('jumps_detected:         ');
while ti < params.tfin
    fprintf('\b\b\b\b\b\b\b%5.01f%%\n',ti/params.ttot*100); 

    % First check if current fstat goes over threshold
    % if not, increment position otherwise proceed with detection
    if Fstats(ti) < params.Fstat_thresh
        ti = ti+1;
    else
        % Once you found the start of the jump find the end by either 
        % finding the point in which Fstat goes below the threshold, or
        % end of timeseries
        tii = 0;
        while Fstats(ti+tii) > params.Fstat_thresh && tii < params.tfin 
            tii = tii+1;
        end
        t_above_thresh = tii;

        % Calculate statistics about each detected jump 
        % columns are: Fstat max, time index of Fstat max, 
        % time index of peak estimated from left edge,
        % time index of peak estimated going backwards from right edge
        peakstats_simple = get_peak_stats_simple(Fstats,ti,ti+tii,params.tfin); 
        Fstatmax = peakstats_simple(1);
        % by convention jump happens 1 sample before it can be measured


        % Set the jump time based off of a particular vaalue in peak stats
        if strcmp(params.jump_type,'inst')
            ti_jump = peakstats_simple(2)+params.Njump_offset;;
        else
            ti_jump = peakstats_simple(4)+params.Njump_offset;;
        end

        % Store relevant infromatino about the jump to return
        % columns are time of detected jump, max of F statistic of jump, 
        % time index of jump, time index crossing above then below threshold
        jumps_detected = [jumps_detected; params.tvect(ti_jump) Fstatmax ti_jump ti ti+tii];
                      
        % start looking for next jump after this one
        ti = ti + tii + 1;
    end
end