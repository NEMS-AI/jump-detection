function Fstat = calc_fstat(fvect_x, fvect_y, params)
% For two datasets, calculate the Fstat between them
% relies on having the nmodes and nmeasure parameter passed through


xbar = mean(fvect_x,2);
ybar = mean(fvect_y,2);
sigma_x = cov(fvect_x');
sigma_y = cov(fvect_y');
sigma_pool = sigma_x/2 + sigma_y/2;
t2 = params.Nmeas/2*((xbar-ybar)'/sigma_pool)*(xbar-ybar);
p_dim = params.nmodes;
Fstat = (2*params.Nmeas-1)/(p_dim*(2*params.Nmeas-2))*t2;

end