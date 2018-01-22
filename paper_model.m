%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 4);
set(0, 'DefaultAxesFontSize', 24);

%% load profiles
% excitation
fiber_profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
fiber_profile_exc = sp_3d_to_2d(fiber_profile_exc); % symmetric, way faster

% emission
fiber_profile_emi = sp_model('sensitivity-profile/fiber-fluor.mat');
fiber_profile_emi = sp_3d_to_2d(fiber_profile_emi); % symmetric, way faster

%% figure 1: fiber distribution

% constants
number_of_fibers = 250;
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

% generate fibers
old_rng = rng; rng(1);
[fibers, fiber_angles] = generate_fibers(number_of_fibers, ...
    'fiber_distribution', distribution, 'angle_distribution', 0, 'volume', volume, ...
    'position', position);
cells = generate_cells('volume', volume);
rng(old_rng);

% calculate center
center = volume .* position;

% calculate cell mixing matrix
fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, cells, ...
    fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
cell_bright = max(fiber_mix, [], 1);
[cell_bright_sort, idx] = sort(cell_bright ./ norm);

% normalization
%norm = max(fiber_profile_emi.volume(:)) * max(fiber_profile_exc.volume(:)); % based on a single fiber maximum potential round-trip brightness
%norm = max(cell_bright); % based on all fibers, brightest cell
norm = max(fiber_profile_emi.volume(:)) * max(fiber_profile_exc.volume(:));

% open window
figure;

% plot cells
threshold = 0.01;
idx = idx(cell_bright_sort >= threshold);
cell_bright_sort = cell_bright_sort(cell_bright_sort >= threshold);
scatter(cells(1, idx) - center(1), cells(2, idx) - center(2), 128, log10(cell_bright_sort), 'filled');
caxis([-2 0]);
colormap('parula');
h = colorbar('Ticks', [-2 -1 0], 'TickLabels', {'10^{-2}', '10^{-1}', '10^0'});
ylabel(h, 'Fluorescence yield [% of max]');
title(sprintf('Fibers: %d; Cells: %d', number_of_fibers, length(cell_bright_sort)));

% plot fibers
hold on;
plot(fibers(1, :) - center(1), fibers(2, :) - center(2), '.', 'MarkerEdgeColor', [1, 0.0784, 0.576], 'MarkerSize', 12);
hold off;
axis square;
l = max(max(abs(fibers(1, :) - center(1))), max(abs(fibers(2, :) - center(2))));
l = 500; % ceil(l / 100) * 100;
xlim([-l l]); xticks([-l 0 l]); xlabel('x [{\mu}]');
ylim([-l l]); yticks([-l 0 l]); ylabel('y [{\mu}]');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig1-distribution.eps'); close;

%% figure 2: sample signals

close all;

% constants
number_of_fibers = 100;
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [25 0 0; 0 25 0; 0 0 5];

old_rng = rng; rng(0);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
    'duration', 200, 'number_of_outputs', number_of_fibers, 'output_noise', @(n) normrnd(0, 0.05, 1, n), ... %
    'params_fibers', {'fiber_distribution', distribution, 'volume', volume, 'position', position}, ...
    'params_cells', {'volume', volume}, 'figures', 2);
rng(old_rng);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig2-signals.eps'); close;
close all;

%% figure 3: auc

close all;

old_rng = rng; rng(0);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, 'duration', 200, 'number_of_outputs', 100, 'auc_threshold', 0.2:0.1:0.8);
rng(old_rng);

h = figure(7);
h.Position = [1000 200 880 812];
r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig4-auc.eps'); close;
close all;

%% figure 4: cell density vs fiber density
explore({'mode', 'profile-rt', 'duration', 200, 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, 'params_fibers', {'angle_distribution', 0, 'fiber_distribution', [100 0 0; 0 100 0; 0 0 0], 'distribution', 'uniform'}}, ...
    'number_of_outputs', 50:50:500, ...
    'param_b_name', param_b_values, 5, true);
