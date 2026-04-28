function plotOverlay()

    % load and normalize experimental
    [exp_data, exp_f1, exp_f2] = experimental(false);
    exp_norm = exp_data / max(abs(exp_data(:)), [], 'omitnan');
    exp_norm(isnan(exp_norm)) = 0;
  
    % sim
    C = setC(zeros(3), [20.0, 29.0], 0); 
    [sim_spec, sim_f1, sim_f2] = simAlanine(C, false);
    sim_spec = imgaussfilt(sim_spec, 1.75);

    figure('Name', 'Exp vs Sim Overlay');
    % plot exp 
    contour(exp_f2, exp_f1, exp_norm, 20, 'r'); hold on;
    % plot sim
    contour(sim_f2, sim_f1, sim_spec, 20, 'k');

    set(gca, 'XDir', 'reverse', 'YDir', 'reverse');
    title('Experimental vs. Simulated');
    xlabel('F2 (ppm)'); ylabel('F1 (ppm)');
    legend('Experimental (Red)', 'Simulated (Black)');
    grid on;

end