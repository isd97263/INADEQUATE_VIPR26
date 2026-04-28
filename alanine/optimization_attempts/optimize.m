function best_J = optimize()

    % load and normalize experimental
    [exp_spec, f1_exp, f2_exp] = experimental(false);
    exp_norm = exp_spec / max(abs(exp_spec(:)), [], 'omitnan');
    exp_norm(isnan(exp_norm)) = 0;

    % hyperparams
    learning_rate_init = 1.0;
    momentum = 0.6;
    max_iterations = 250;
    h_eff = 0.3;      % +/- .3 hz  

    tolerance = 1e-7;
    patience = 10;
    num_restarts = 5;
    best_overall_cost = inf;
    
    lb = [10.0, 10.0, 0.0]; % lower bounds for coupling [J12, J23, J13]
    ub = [30.0, 40.0, 2.0]; % upper bounds
    
    log_file = 'j_all_progress.txt';
    fid = fopen(log_file, 'w');
    fprintf(fid, 'Global 3-Parameter Optimization Started: %s\n', datestr(now));
    fclose(fid);

    for restart = 1:num_restarts

        % Randomly initialize all 3 parameters 
        current_J = lb + (ub - lb) .* rand(1, 3); 
        velocity = [0, 0, 0];
        
        restart_best_cost = inf;
        restart_best_J = current_J;

        prev_cost = inf;
        stall_counter = 0;

        fid = fopen(log_file, 'a');
        fprintf(fid, '\n--- RESTART %d/%d: Initial [J12, J23, J13] = [%.3f, %.3f, %.3f] Hz ---\n', ...
            restart, num_restarts, current_J(1), current_J(2), current_J(3));
        fclose(fid);

        for iter = 1:max_iterations

            calc_cost = @(J_vec) compute_cost(J_vec(1), J_vec(2), J_vec(3), exp_norm, f1_exp, f2_exp);

            current_cost = calc_cost(current_J);

            if current_cost < restart_best_cost
                restart_best_cost = current_cost;
                restart_best_J = current_J;
            end

            grad = zeros(1, 3);
            for p = 1:3
                % Perturb current parameter
                J_plus = current_J; J_plus(p) = J_plus(p) + h_eff;
                J_minus = current_J; J_minus(p) = J_minus(p) - h_eff;
                
                cost_plus = calc_cost(J_plus);
                cost_minus = calc_cost(J_minus);
                
                grad(p) = (cost_plus - cost_minus) / (2 * h_eff);
            end

            current_lr = learning_rate_init * (1 / (1 + 0.01 * iter));
            velocity = momentum * velocity - current_lr * grad;
            current_J = current_J + velocity;

            % Enforce boundaries on all 3 parameters 
            current_J = max(lb, min(ub, current_J));

            fid = fopen(log_file, 'a');
            fprintf(fid, 'Iter %3d | J12: %6.5f | J23: %6.5f | J13: %6.5f | Cost: %.8f \n', ...
                iter, current_J(1), current_J(2), current_J(3), current_cost);
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
            best_overall_J = restart_best_J;
        end
    end

    best_J = best_overall_J;

    fid = fopen(log_file, 'a');
    fprintf(fid, '\n========================================\n');
    fprintf(fid, 'OPTIMIZATION COMPLETE.\n');
    fprintf(fid, 'Global Best J12: %.4f Hz\n', best_J(1));
    fprintf(fid, 'Global Best J23: %.4f Hz\n', best_J(2));
    fprintf(fid, 'Global Best J13: %.4f Hz\n', best_J(3));
    fprintf(fid, 'Final Cost: %.5f\n', best_overall_cost);
    fprintf(fid, '========================================\n');
    fclose(fid);
    
    save('j_all_simple_result.mat', 'best_J', 'best_overall_cost');
end


function cost = compute_cost(j12, j23, j13, exp_norm, f1_exp, f2_exp)
    try
        C = setC(zeros(3), [j12, j23], j13);
        [sim_norm, sim_f1, sim_f2] = simAlanine(C, false);

        [X, Y] = meshgrid(flip(f2_exp), flip(f1_exp));

        sim_aligned_asc = interp2(flip(sim_f2), flip(sim_f1), flipud(fliplr(sim_norm)), X, Y, 'spline', 0);
        sim_aligned = flipud(fliplr(sim_aligned_asc));
        
        
        mse = mean((sim_aligned(:) - exp_norm(:)).^2) * 1000000;

        cost = mse;
        
    catch
        cost = 1e6;
    end
end