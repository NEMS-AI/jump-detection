function peak_stats = get_peak_stats(Fstats,tvect,ti1,ti2,tfin)
%get_peak_stats summarize various statistics related to Fstat peak.
% Fstats input is the variable over entire dataset
% tvect is similarly time vector over entire dataset
% ti1: time index of Fstats crossing above threshold
% ti2: time index of Fstats crossing back below threshold
% tfin: time index of end of dataset (TODO: not needed?)

    ti_start = ti1;
    ti_end = ti2;

    Fstatmax = max(Fstats(ti_start:ti_end));
    ti_Fstatmax = find(Fstats(ti_start:ti_end)==Fstatmax)+ti_start-1;
    
    % find time and F stat vector from start of first jump window to
    % time of jump
    tvect_i = tvect(ti_start:ti_end);
    Fstats_i = log10(Fstats(ti_start:ti_end)); % take log
%     alternative normalization
%     Fstats_i = Fstats(ti_start:ti_end)/max(Fstats(ti_start:ti_end));
    
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
    t_above_thresh = 0;
    
    % normalize Fstats and collect statistics
    Fstats_i = max(Fstats_i,Fstats_i(end));
    Fstats_i = Fstats_i-Fstats_i(end);
    Fint_0 = trapz(tvect_i,Fstats_i);
    t_mean = trapz(tvect_i,tvect_i.*Fstats_i/Fint_0);
    t_mean_rel = t_mean - tvect(ti_end);
    F_var = trapz(tvect_i,(tvect_i-t_mean).^2.*Fstats_i/Fint_0);  
    F_skew = trapz(tvect_i,(tvect_i-t_mean).^3.*Fstats_i/Fint_0)/F_var^1.5;
    F_kurt = trapz(tvect_i,(tvect_i-t_mean).^4.*Fstats_i/Fint_0)/F_var^2;
    F_var2 = trapz(tvect_i,(tvect_i-t_mean).^2.*Fstats_i/max(Fstats_i));
    F_kurt2 = trapz(tvect_i,(tvect_i-t_mean).^4.*Fstats_i/max(Fstats_i));

    peak_stats = [Fstatmax t_mean_rel sqrt(F_var) F_skew F_kurt sqrt(F_var2) F_kurt2 tvect(t_fwhm+1)];
end

