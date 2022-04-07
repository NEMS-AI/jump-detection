function peak_stats = get_peak_stats(Fstats,tvect,ti1,ti2,ti3,ti4,ti5,tfin)
%get_peak_stats summarize various statistics related to Fstat peak.
% Fstats input is the variable over entire dataset
% tvect is similarly time vector over entire dataset
% ti1: time index of beginning of pre-jump measurement window
% ti2: time index of end of pre-jump measurement window (unused right now)
% ti3: time index of end of post-jump measurement window
% ti4: time index of Fstats crossing above threshold
% ti5: time index of Fstats crossing back below threshold
% tfin: time index of end of dataset (TODO: not needed?)
% Fstat is above threshold between ti4 and ti5.

    Fstatmax = max(Fstats(ti1:ti5));
    ti_Fstatmax = find(Fstats(ti1:ti5)==Fstatmax)+ti1-1;
    
    % find time and F stat vector from start of first jump window to
    % time of jump
    tvect_i = tvect(ti1:ti5);
%     Fstats_i = Fstats(ti1:ti3)/max(Fstats(ti1:ti3)); % normalize
    Fstats_i = log10(Fstats(ti1:ti5)); % take log
    
    peaks_right = 0; % flag if F stat goes significantly above threshold in post-jump window; not used
    tvect_right = tvect(ti5+1:ti3); % time and F stats vector after jump time
    Fstats_right = log10(Fstats(ti5+1:ti3));
    if sum(Fstats_right > Fstats_i(end)*1.25) > 0, peaks_right = 1; end
    
    %  find FWHM 
    peak_left_i = ti_Fstatmax;
    while Fstats(peak_left_i) > Fstatmax/2 && peak_left_i > 1
        peak_left_i = peak_left_i - 1;
    end
    peak_right_i = ti_Fstatmax;
    while Fstats(peak_right_i) > Fstatmax/2 && peak_right_i < tfin
        peak_right_i = peak_right_i + 1;
    end
    t_fwhm = peak_right_i - peak_left_i;
    
    % find time above threshold for primary peak
    peak_left_i = ti_Fstatmax;
    while Fstats(peak_left_i) > Fstats(ti5) && peak_left_i > ti1
        peak_left_i = peak_left_i - 1;
    end
    t_above_thresh = ti5 - peak_left_i;
    
    % normalize Fstats and collect statistics
    Fstats_i = max(Fstats_i,Fstats_i(end));
    Fstats_i = Fstats_i-Fstats_i(end);
    Fint_0 = trapz(tvect_i,Fstats_i);
    t_mean = trapz(tvect_i,tvect_i.*Fstats_i/Fint_0);
    t_mean_rel = t_mean - tvect(ti5);
    F_var = trapz(tvect_i,(tvect_i-t_mean).^2.*Fstats_i/Fint_0);
%     F_skew = trapz(tvect_i,(tvect_i-t_mean).^3.*Fstats_i/Fint_0)/F_var^1.5;
%     F_kurt = trapz(tvect_i,(tvect_i-t_mean).^4.*Fstats_i/Fint_0)/F_var^2;

    t_bif = 0; % time spent above bifurcation point (not used)
    peak_stats = [Fstatmax t_mean_rel tvect(t_fwhm+1) tvect(t_bif+1) tvect(t_above_thresh+1) sqrt(F_var) peaks_right];
end

