function compare_mse() 
    
    % load and normalize experimental
    [exp_data, f1_exp, f2_exp] = experimental(false);
    exp_norm = exp_data / max(abs(exp_data(:)), [], 'omitnan');
    exp_norm(isnan(exp_norm)) = 0;

    C = setC(zeros(3), [35.31/2, 54.979/2], 0);
    [sim_actual, sim_f1, sim_f2] = simAlanine(C, false);

    % match native physical axis grid of experimental
    % original matrix elements do not map onto each other even if
    % visually do
    [X, Y] = meshgrid(flip(f2_exp), flip(f1_exp));

    sim_aligned_asc = interp2(flip(sim_f2), flip(sim_f1), flipud(fliplr(sim_actual)), X, Y, 'spline', 0);
    sim_aligned = flipud(fliplr(sim_aligned_asc));

    mse_actual = mean((sim_aligned(:) - exp_norm(:)).^2) * 100000000;

    C = setC(zeros(3), [17.655, 20.4284], 0);
    [sim_opt, sim_f1, sim_f2] = simAlanine(C, false);

    % match native physical axis grid of experimental
    % original matrix elements do not map onto each other even if
    % visually do
    [X, Y] = meshgrid(flip(f2_exp), flip(f1_exp));

    sim_aligned_asc = interp2(flip(sim_f2), flip(sim_f1), flipud(fliplr(sim_opt)), X, Y, 'spline', 0);
    sim_aligned = flipud(fliplr(sim_aligned_asc));

    mse_opt = mean((sim_aligned(:) - exp_norm(:)).^2) * 100000000;

    fprintf('ACTUALS | MSE: %.10f\n', mse_actual);
    fprintf('OPTIMIZED | MSE: %.10f\n', mse_opt);
    disp('=============================================');

end