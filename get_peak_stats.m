function peak_stats = get_peak_stats(Fstats,tvect,ti1,ti2,tfin)
%get_peak_stats summarize various statistics related to Fstat peak.
%Assume Fstat is above threshold between ti1 and ti2

    ti = ti1; tii = ti2-ti1;
    t_above_thresh = tii;
    Fstatmax = max(Fstats(ti:ti+tii));
    ti_Fstatmax = find(Fstats(ti:ti+tii)==Fstatmax);
    
    %  find approx FWHM 
    peak_left_i = ti+ti_Fstatmax;
    while Fstats(peak_left_i) > Fstatmax/2 && peak_left_i > 1
        peak_left_i = peak_left_i - 1;
    end
    peak_right_i = ti+ti_Fstatmax;
    while Fstats(peak_right_i) > Fstatmax/2 && peak_right_i < tfin
        peak_right_i = peak_right_i + 1;
    end
    t_fwhm = peak_right_i - peak_left_i;

    % find mean, std, skew.
    tvect_i = tvect(ti:ti+tii);
    Fstats_i = Fstats(ti:ti+tii);
    Fint_0 = trapz(tvect_i,Fstats_i);
    t_mean = trapz(tvect_i,tvect_i.*Fstats_i/Fint_0);
    F_var = trapz(tvect_i,(tvect_i-t_mean).^2.*Fstats_i/Fint_0);
    F_skew = trapz(tvect_i,(tvect_i-t_mean).^3.*Fstats_i/Fint_0)/F_var^1.5;
    F_kurt = trapz(tvect_i,(tvect_i-t_mean).^4.*Fstats_i/Fint_0)/F_var^2;
    
    peak_stats = [Fstatmax Fint_0 t_mean tvect(t_above_thresh+1) tvect(t_fwhm+1) sqrt(F_var) F_skew F_kurt];
end

