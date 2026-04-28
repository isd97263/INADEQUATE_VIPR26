function best_j23 = optimize_j23_mse(best_j12)

    % load and normalize experimental
    [exp_spec, f1_exp, f2_exp] = experimental(false);
    exp_norm = exp_spec / max(abs(exp_spec(:)), [], 'omitnan');
    exp_norm(isnan(exp_norm)) = 0;

    % hyperparams
    learning_rate_init = 1.0;
    momentum = 0.6;
    max_iterations = 250;
    h_eff = 0.2;      % +/- .2 hz  

    tolerance = 1e-5;
    patience = 10;

    num_restarts = 5;
    best_overall_cost = inf;

    % fixed parameters for alanine
    J12_fixed = best_j12;
    J13_fixed = 0.0;
    
    log_file = 'j23_simple_progress.txt';
    fid = fopen(log_file, 'w');
    fprintf(fid, 'J23 Optimization Started: %s\n', datestr(now));
    fprintf(fid, '========================================\n');
    fclose(fid);

    for restart = 1:num_restarts

        % randomly initialize J23
        current_j23 = 10 + (30 - 10) * rand();
        velocity = 0;
        restart_best_cost = inf;

        prev_cost = inf;
        stall_counter = 0;

        fid = fopen(log_file, 'a');
        fprintf(fid, '\n--- RESTART %d/%d: Initial J23 = %.3f Hz ---\n', restart, num_restarts, current_j23);
        fclose(fid);

        for iter = 1:max_iterations

            calc_cost = @(j) compute_cost(J12_fixed, j, J13_fixed, exp_norm, f1_exp, f2_exp);

            current_cost = calc_cost(current_j23);

            if current_cost < restart_best_cost
                restart_best_cost = current_cost;
            end

            cost_plus = calc_cost(current_j23 + h_eff);
            cost_minus = calc_cost(current_j23 - h_eff);
            grad = (cost_plus - cost_minus) / (2 * h_eff);

            current_lr = learning_rate_init * (1 / (1 + 0.01 * iter));
            velocity = momentum * velocity - current_lr * grad;
            current_j23 = current_j23 + velocity;

            current_j23 = max(10, min(30, current_j23));

            fid = fopen(log_file, 'a');
            fprintf(fid, 'Iter %3d | J23: %6.3f Hz | Cost: %.5f \n', ...
                iter, current_j23, current_cost);
            fclose(fid);

             cost_diff = abs(prev_cost - current_cost);
            if cost_diff < tolerance
                stall_counter = stall_counter + 1;
            else 
                stall_counter = 0;
            end
            prev_cost = current_cost;

            if stall_counter >= patience
                fid = fopen(log_file, 'a');
                fprintf(fid, 'Early stopping at iteration %d \n', iter);
                fclose(fid);
                break;
            end

        end

        if restart_best_cost < best_overall_cost
            best_overall_cost = restart_best_cost;
            best_overall_j23 = current_j23;
        end
    end

    best_j23 = best_overall_j23;

    fid = fopen(log_file, 'a');
    fprintf(fid, '\n========================================\n');
    fprintf(fid, 'OPTIMIZATION COMPLETE.\n');
    fprintf(fid, 'Global Best J23: %.4f Hz\n', best_j23);
    fprintf(fid, 'Final Cost: %.5f\n', best_overall_cost);
    fprintf(fid, '========================================\n');
    fclose(fid);
    
    save('j23_simple_result.mat', 'best_j23', 'best_overall_cost');
end


function cost = compute_cost(j12, j23, j13, exp_norm, f1_exp, f2_exp)
    try

        C = setC(zeros(3), [j12, j23], j13);
        [sim_norm, sim_f1, sim_f2] = simAlanine(C, false);

        % match native physical axis grid of experimental
        % original matrix elements do not map onto each other even if
        % visually do
        [X, Y] = meshgrid(flip(f2_exp), flip(f1_exp));

        sim_aligned_asc = interp2(flip(sim_f2), flip(sim_f1), flipud(fliplr(sim_norm)), X, Y, 'spline', 0);
        sim_aligned = flipud(fliplr(sim_aligned_asc));

        mse = mean((sim_aligned(:) - exp_norm(:)).^2) * 1000;
        cost = mse
        
    catch
        cost = 1e6;
    end
end