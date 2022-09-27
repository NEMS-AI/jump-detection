function Fstats = calculate_fstats(params)
% Calculate Fstats Based of input parameters, of a time series, calculate
% Fstats


Fstats = zeros(params.tfin,1);
fprintf('Fstats:         ');
for ti=1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    :params.tfin-params.tlag_buffer % allow for time lag in either direction
   
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',ti/params.ttot*100); 
    
    xi_range = ti:ti+params.Nmeas-1;  
    yi_range = ti+params.Nmeas+params.Njump:ti+2*params.Nmeas+params.Njump-1;

    fvect_xi = params.fvect(:,xi_range);
    fvect_yi = params.fvect(:,yi_range);

    Fstat = calc_fstat(fvect_xi,fvect_yi, params);
    Fstats(ti+params.Nmeas) = Fstat;
    
end

% plot Fstats
figure;
plot(params.tvect(1:length(Fstats)),Fstats);
xlabel('Time (s)')
ylabel('F statistic');

end
