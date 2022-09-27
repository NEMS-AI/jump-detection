function [rel_jump, rel_err, err_rho, xbar, ybar] = calculate_jump(params, fvect_x, fvect_y)
    xbar = median(fvect_x,2);
    ybar = median(fvect_y,2);
    sigma_x = cov(fvect_x');
    sigma_y = cov(fvect_y');
    sigma_pool = sigma_x/2 + sigma_y/2;
    rel_jump = ybar./xbar-1;

    % TODO Redo this for-loop as matrix operations
    rel_err = [];
    for i = 1:params.nmodes
    rel_err = [rel_err sqrt(sigma_pool(i,i))];
    end
    rel_err = rel_err./params.fvect(:,1)';
    
    % TODO Fix this for nmode 3
    err_rho = sigma_pool(1,2)/(sqrt(sigma_pool(1,1))*sqrt(sigma_pool(2,2)));
end