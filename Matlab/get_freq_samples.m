%% get frequency samples (fvect_xi, fvect_yi) for given time and time lag per second
tlag_m2 = (tvect(ti)+tstart)*tlag_mode2_per_second;
Nlag_m2 = floor(tlag_m2/tsample);
fvect_xi = [fvect(1,xi_range); fvect(2,xi_range+Nlag_m2)];
fvect_yi = [fvect(1,yi_range); fvect(2,yi_range+Nlag_m2)];

if nmodes == 3
    tlag_m3 = (tvect(ti)+tstart)*tlag_mode3_per_second;
    Nlag_m3 = floor(tlag_m3/tsample);
    fvect_xi = [fvect_xi; fvect(3,xi_range+Nlag_m3)];
    fvect_yi = [fvect_yi; fvect(3,yi_range+Nlag_m3)];
end