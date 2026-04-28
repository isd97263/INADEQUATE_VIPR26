%{
Returns the experimental spectrum that is masked to only show peak centers
and referenced to the top right peak at (18.62, 71.48)
%}

function [referenced_spec, f1_axis, f2_axis] = experimental(plot)

    experimental_spec = Load2D('/home/nmrbox/0054/idelgado/2025/experimental/alanine');

    masked_spec = experimental_spec.real;

    f1_axis = experimental_spec.ppm2;
    f2_axis = experimental_spec.ppm1;

    peak_centers = [
        175.75, 126.2;
        50.45, 126.2;
        50.45, 91.8;
        50.45, -33.5;
        16.06, -33.5;
    ];

    f2_width = 5;
    f1_width = 10;

    final_mask = false(size(masked_spec));

    for i = 1:size(peak_centers, 1)
        center_f2 = peak_centers(i, 1);
        center_f1 = peak_centers(i, 2);

        % Find indices within ranges 
        f1_range = abs(f1_axis - center_f1) <= f1_width;
        f2_range = abs(f2_axis - center_f2) <= f2_width;

        % Create the 2D box 
        peak_box = (f1_range(:) & f2_range(:).');

        final_mask = final_mask | peak_box;
    end

    masked_spec(~final_mask) = NaN;

    % target = top right peak
    current_f2 = 16.065;
    target_f2 = 18.62;

    current_f1 = -33.52;
    target_f1 = 71.48;

    ppm_per_point_f2 = abs(f2_axis(1) - f2_axis(2));
    ppm_per_point_f1 = abs(f1_axis(1) - f1_axis(2));
    shift_f2_points = round((current_f2 - target_f2) / ppm_per_point_f2);
    shift_f1_points = round((current_f1 - target_f1) / ppm_per_point_f1);
    referenced_spec = circshift(masked_spec, [shift_f1_points, shift_f2_points]);

    if plot
        max_int = max(referenced_spec(:), [], 'omitnan');
        my_levels = linspace(0.05*max_int, max_int, 10);
    
        [M, c] = contour(f2_axis, f1_axis, referenced_spec, my_levels);
        c.LineColor = 'r';
        set(gca, 'XDir', 'reverse', 'YDir', 'reverse');
        xlabel('F2: 13C chemical shift / ppm');
        ylabel('F1: DQ dimension / ppm');
        ylim([0 250]);
        hold off;
    end

end

