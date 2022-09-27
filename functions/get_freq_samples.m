%% get frequency samples (fvect_xi, fvect_yi) for given time and time lag per second
function [fvect_xi,fvect_yi ] = get_freq_samples(params, ti)

tlag_m2 = (params.tvect(ti)+params.tstart)*params.tlag_mode2_per_second;
Nlag_m2 = floor(tlag_m2/params.tsample);
fvect_xi = [params.fvect(1,xi_range); params.fvect(2,xi_range+Nlag_m2)];
fvect_yi = [params.fvect(1,yi_range); params.fvect(2,yi_range+params.Nlag_m2)];

end