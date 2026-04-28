function combined_spectrum = simGlucose(alpha_coupling, beta_coupling)

forms = {'beta', 'alpha'};

zeeman_sets = {
    [63.28, 74.11, 72.32, 75.48, 74.16, 94.79],  % a-glucose
    [63.43, 78.63, 72.30, 78.44, 76.82, 98.61]     % b-glucose
};

parameters.spins={'13C'};
parameters.offset = 22132.85; %100 ppm
parameters.sweep = [13888.889 27668.273]; %[47169.811, 94074.009]; %to match frac 55 --> [170, 340]
parameters.npoints = [2048 4096];
parameters.zerofill = [4096 8192];
parameters.axis_units = 'ppm';

% Parameters for plotting (both need DQ and SQ offsets)
parameters_plot = parameters;
parameters_plot.offset = [22132.85-1512.8182 - 13.5575 22132.85-735.7547 - 20];

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Allocate storage for both spectra
spectrum_all = cell(1, 2);

parameters.J = 50;

% System setup
sys.isotopes = {'13C','13C','13C','13C','13C','13C'};
sys.magnet = 25.83793;

for f = 1:length(forms) % alpha, betac

    inter.zeeman.scalar = num2cell(zeeman_sets{f});

    if f==1
       inter.coupling.scalar = num2cell(alpha_coupling);
    else
       inter.coupling.scalar = num2cell(beta_coupling);
    end

    inter.temperature = 300;

    % Spinach setup
    spin_system = create(sys, inter);
    spin_system.sys.output = 'hush';
    spin_system = basis(spin_system, bas);
    
    subsystems = {spin_system};
    disp('Using full spin system.');

    % Output prealloc
    spectrum = zeros(parameters.zerofill(2), parameters.zerofill(1));

    % Loop over all isotopomers
    for n = 1:numel(subsystems)
        subsystem = basis(subsystems{n}, bas);
        fid = liquid(subsystem, @inadequate_2d, parameters, 'nmr');
        fid.cos = apodisation(spin_system, fid.cos, {{'cos'}, {'cos'}});
        fid.sin = apodisation(spin_system, fid.sin, {{'cos'}, {'cos'}});
        spec_cos = fftshift(fft(fid.cos, parameters.zerofill(2), 1), 1);
        spec_sin = fftshift(fft(fid.sin, parameters.zerofill(2), 1), 1);
        spec_states = real(spec_cos) + 1i * real(spec_sin);
        spectrum = spectrum + fftshift(fft(spec_states, parameters.zerofill(1), 2), 2);
    end


    % Normalize and store
    spectrum_all{f} = abs(spectrum) / max(abs(spectrum(:)));

    % Individual plot
    % fig = figure('Visible','off'); 
    % scale_figure([1.5 1.5]);
    % plot_2d(spin_system, real(spectrum), parameters_plot, 20, [0.01 0.5 0.01 0.5], 2, 256, 6, 'positive');
    % title(sprintf('INADEQUATE: %s-glucose', forms{f}));
    % kylabel('F1: DQ dimension / ppm');

end % forms

% alpha and beta network overlay
combined_spectrum = spectrum_all{1} + spectrum_all{2};
combined_spectrum = combined_spectrum / max(combined_spectrum(:));  % Normalize

% parameters_plot.offset = [20751.205 20751.205]; %[22132.85 - 790, 22132.85 - 790]; %[35969.862 29079.618]; %[29079.618 35969.862]; % [22132.85, 22132.85]; 
% parameters_plot.sweep = [13888.889 27668.273];
% parameters_plot.npoints = parameters.npoints;


% Overlay plot using plot_2d
fig = figure('Visible','off');
% scale_figure([1.5 1.5]);
[axis_f1, axis_f2, ~] = plot_2d(spin_system, real(combined_spectrum), parameters_plot, ...
     20, [0.01 0.5 0.01 0.5], 2, 256, 6, 'positive');
% title('INADEQUATE Overlay: a- and b-glucose');
% kylabel('F1: DQ dimension / ppm');
% set(gca, 'XDir', 'reverse');
% set(gca, 'YDir', 'reverse');
% xlim([55 105]);
% ylim([30 130]);
close(fig);

figure();
combined_spectrum = combined_spectrum';
[M, c] = contour(axis_f2, axis_f1, combined_spectrum, 20);
c.LineColor = 'b';
set(gca, 'XDir', 'reverse', 'YDir', 'reverse');
xlabel('F2: 13C chemical shift / ppm');
ylabel('F1: DQ dimension / ppm');
xlim([55 105]);
ylim([30 130]);
hold on;
end

