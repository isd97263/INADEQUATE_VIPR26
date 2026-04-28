% match doublets found in small region using csv of peak peaked data

function doublet_data = find_doublet_coupling(f1_min, f1_max, f2_min, f2_max)
    file = 'Frac55INADPeaks_unmerged.csv';
    
    data = readtable(file);
    f1_ppm = data.f1_ppm_;
    f2_ppm = data.f2_ppm_;
    freq_MHz = 276.68;
    f2_Hz = f2_ppm * freq_MHz;
    
    target_region = (f1_ppm >= f1_min & f1_ppm <= f1_max) & ...
                (f2_ppm >= f2_min & f2_ppm <= f2_max);
    
    f1_filtered = f1_ppm(target_region);
    f2_filtered_Hz = f2_Hz(target_region);
    
    J_expected_min = 20;
    J_expected_max = 60;
    
    f1_tolerance = 0.05;

    doublet_data = [];
    processed_f1_groups = []; % f1 coords the loop has already encountered
    
    for i = 1:length(f1_filtered)
      current_f1 = f1_filtered(i); % specific peak

        if any(abs(processed_f1_groups - current_f1) <= f1_tolerance)
            continue;
        end

        processed_f1_groups(end+1) = current_f1;
        
        % all peaks within current peak and tolerance range
        slice_idx = abs(f1_filtered - current_f1) <= f1_tolerance;
        f2_slice_Hz = sort(f2_filtered_Hz(slice_idx));
        
        if length(f2_slice_Hz) == 2
            f2_diff = abs(f2_slice_Hz(2) - f2_slice_Hz(1));
            
            if f2_diff >= J_expected_min && f2_diff <= J_expected_max
                doublet_data(end+1, :) = [f2_diff, current_f1, f2_slice_Hz(1), f2_slice_Hz(2)];
            end
        end
    end
    
    if isempty(doublet_data)
        disp('No valid doublets found in the specified window.');
    else
        num_doublets = size(doublet_data, 1);
        fprintf('Found %d valid doublet(s):\n', num_doublets);
        for k = 1:num_doublets
            fprintf('  Doublet %d\n', k);
            fprintf('  J-coupling:  %.2f Hz\n', doublet_data(k, 1));
            fprintf('  F1 (Shared): %.2f ppm\n', doublet_data(k, 2));
            fprintf('  F2 (Peak 1): %.2f ppm\n', doublet_data(k, 3)/freq_MHz);
            fprintf('  F2 (Peak 2): %.2f ppm\n\n', doublet_data(k, 4)/freq_MHz);
        end
    end
end