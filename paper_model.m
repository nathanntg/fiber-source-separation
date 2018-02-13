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
number_of_fibers = [250 500];
cell_density = [0.000275 0.00078];
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

% open window
h = figure;
h.Position = h.Position .* [1 1 2 2];

% seed random number generator
old_rng = rng; rng(1);

for i = 1:length(number_of_fibers)
    for j = 1:length(cell_density)
        subplot(length(cell_density), length(number_of_fibers), (j - 1) * length(cell_density) + i);

        % generate fibers
        [fibers, fiber_angles] = generate_fibers(number_of_fibers(i), ...
            'fiber_distribution', distribution, 'angle_distribution', 0, 'volume', volume, ...
            'position', position);
        cells = generate_cells('volume', volume, 'cell_density', cell_density(j));

        % calculate center
        center = volume .* position;

        % normalization
        %norm = max(fiber_profile_emi.volume(:)) * max(fiber_profile_exc.volume(:)); % based on a single fiber maximum potential round-trip brightness
        %norm = max(cell_bright); % based on all fibers, brightest cell
        norm = max(fiber_profile_emi.volume(:)) * max(fiber_profile_exc.volume(:));

        % calculate cell mixing matrix
        fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, cells, ...
            fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
        cell_bright = max(fiber_mix, [], 1);
        [cell_bright_sort, idx] = sort(cell_bright ./ norm);

        % plot cells
        threshold = 0.01;
        idx = idx(cell_bright_sort >= threshold);
        cell_bright_sort = cell_bright_sort(cell_bright_sort >= threshold);
        scatter(cells(1, idx) - center(1), cells(2, idx) - center(2), 128, log10(cell_bright_sort), 'filled');
        caxis([-2 0]);
        colormap('parula');
        %title(sprintf('Fibers: %d; Neurons: %dk per mm^3', number_of_fibers(i), cell_density(j) * 1000^3 / 1000));

        % plot fibers
        hold on;
        plot(fibers(1, :) - center(1), fibers(2, :) - center(2), '.', 'MarkerEdgeColor', [1, 0.0784, 0.576], 'MarkerSize', 14);
        hold off;
        axis square;
        l = max(max(abs(fibers(1, :) - center(1))), max(abs(fibers(2, :) - center(2))));
        l = 520; % ceil(l / 100) * 100;
        xlim([-l l]); xticks([-500 0 500]); xlabel('x [{\mu}]');
        ylim([-l l]); yticks([-500 0 500]); ylabel('y [{\mu}]');
        
        % number of cells
        h = text(.98 * l, .98 * -l, sprintf('Neurons: %d', length(cell_bright_sort)), 'FontSize', 22, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

% colorbar position
spp1 = get(subplot(length(cell_density), length(number_of_fibers), 1), 'Position');
spp2 = get(subplot(length(cell_density), length(number_of_fibers), 2), 'Position');
spp3 = get(subplot(length(cell_density), length(number_of_fibers), 3), 'Position');
spp4 = get(subplot(length(cell_density), length(number_of_fibers), 4), 'Position');
h = colorbar('Ticks', [-2 -1 0], 'TickLabels', {'10^{-2}', '10^{-1}', '10^0'}, ...
    'Position', [spp2(1) + spp2(3), spp4(2), 0.02, spp2(2) + spp2(4) - spp4(2)]);
ylabel(h, 'Fluorescence yield [frac. of max]');

% titles
ax = axes('Position', [0 0 1 1], 'Visible', 'off');
text(ax, spp1(1) + spp1(3) / 2, spp1(2) + spp1(4) + 0.02, sprintf('Fibers: %d', number_of_fibers(1)), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized');
text(ax, spp2(1) + spp2(3) / 2, spp2(2) + spp2(4) + 0.02, sprintf('Fibers: %d', number_of_fibers(2)), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized');
text(ax, spp1(1) - 0.07, spp1(2) + spp1(4) / 2, sprintf('Density: %dk per mm^3', cell_density(1) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized', 'Rotation', 90);
text(ax, spp3(1) - 0.07, spp3(2) + spp3(4) / 2, sprintf('Density: %dk per mm^3', cell_density(2) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized', 'Rotation', 90);

rng(old_rng);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig1-distribution.eps'); close;

%% figure 2: number of cells seen

% constants
number_of_fibers = [50 100 150 200 500 750 1000 1500 2000];
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

iters = 5;
thresholds = [0.01 0.02 0.05 0.1];

result_single = zeros(length(number_of_fibers), length(thresholds), iters);
result_multi = zeros(length(number_of_fibers), length(thresholds), iters);
result_per_fiber = zeros(length(number_of_fibers), length(thresholds), iters);

% generate fibers
old_rng = rng; rng(0);
cells = generate_cells('volume', volume);

for i = 1:length(number_of_fibers)
    for j = 1:iters
        % generate fibers
        [fibers, fiber_angles] = generate_fibers(number_of_fibers(i), ...
            'fiber_distribution', distribution, 'volume', volume, ...
            'position', position);
        
        fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, cells, ...
            fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
        cell_bright = max(fiber_mix, [], 1);
        
        for k = 1:length(thresholds)
            % swap threshold and iter index
            result_single(i, k, j) = sum(cell_bright >= thresholds(k));
            result_multi(i, k, j) = sum(sum(fiber_mix >= thresholds(k), 1) >= 2);
            result_per_fiber(i, k, j) = mean(sum(fiber_mix >= thresholds(k), 2));
        end
    end
end

rng(old_rng);

%%

% plot it
h = figure;
h.Position(3) = 2 * h.Position(3);

subplot(1, 2, 1);

mn = mean(result_single, 3);
% st = std(result_single, 0, 3);
set(gca, 'ColorOrderIndex', 1);
h = loglog(number_of_fibers, mn);
% hold on;
% plot(number_of_fibers, mn + st, 'LineWidth', 1);
% plot(number_of_fibers, mn - st, 'LineWidth', 1);
% hold off;
title('Neurons');
xlabel('Number of fibers');
ylabel('Neurons above threshold');

hold on;
mn = mean(result_multi, 3);
% st = std(result_single, 0, 3);
set(gca, 'ColorOrderIndex', 1);
loglog(number_of_fibers, mn, ':');
hold off;

% adjust dimensions
ylim([10 14000]);
xlim(number_of_fibers([1 end]));

hold on;
mn = mean(result_multi, 3);
% st = std(result_single, 0, 3);
set(gca, 'ColorOrderIndex', 1);
plot(number_of_fibers, mn, ':');
hold off;

l = cell(1, length(thresholds));
for i = 1:length(thresholds)
    l{i} = sprintf('\\geq %.1f%%', thresholds(i)*100);
end
legend(h, l{:}, 'Location', 'NorthWest');

subplot(1, 2, 2);

mn = mean(result_per_fiber, 3);
set(gca, 'ColorOrderIndex', 1);
loglog(number_of_fibers, mn);

xlim(number_of_fibers([1 end]));
ylim([0 130]);

xlabel('Number of fibers');
ylabel('Neurons above threshold');
title('Neurons per fiber');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig2-cells.eps'); close;

%% figure 3: neuron strengths per fiber

% constants
number_of_fibers = [50 250 500 750];
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

top_neurons = 10; % top most visible fibers
iters = 5; % number of iterations
norm = true; % normalize by brigtest fiber

result = cell(1, length(number_of_fibers));

% generate fibers
old_rng = rng; rng(0);
cells = generate_cells('volume', volume);

for i = 1:length(number_of_fibers)
    result{i} = zeros(number_of_fibers(i), top_neurons, iters);
    for j = 1:iters
        % generate fibers
        [fibers, fiber_angles] = generate_fibers(number_of_fibers(i), ...
            'fiber_distribution', distribution, 'volume', volume, ...
            'position', position);
        
        fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, cells, ...
            fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
        
        % visbility of top neurons
        m = sort(fiber_mix, 2, 'descend');
        m = m(:, 1:top_neurons);
        if norm
            m = bsxfun(@rdivide, m, m(:, 1));
        else
            m = m ./ max(fiber_profile_emi.volume(:)) * max(fiber_profile_exc.volume(:));
        end
        
        result{i}(:, :, j) = m;
    end
end

rng(old_rng);

h = figure;
h.Position = h.Position .* [1 1 2 2];

for i = 1:length(number_of_fibers)
    subplot(2, 2, i);
    boxplot(reshape(permute(result{i}, [2 1 3]), top_neurons, [])' .* 100);
    xlabel('Neurons in descending order of fluence');
    if norm
        ylabel('% of brightest');
    else
        ylabel('% of max');
    end
    title(sprintf('%d fibers', number_of_fibers(i)));
end

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig3-neurons.eps'); close;

%% figure 4: number of fibers per cell

% constants
number_of_fibers = [50 250 500 750];
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

top_neurons = 50; % top brightest neurons
top_fibers = 10; % top most visible fibers
iters = 5; % number of iterations
norm = true; % normalize by brigtest fiber

result = zeros(length(number_of_fibers), top_fibers, top_neurons, iters);

% generate fibers
old_rng = rng; rng(0);
cells = generate_cells('volume', volume);

for i = 1:length(number_of_fibers)
    for j = 1:iters
        % generate fibers
        [fibers, fiber_angles] = generate_fibers(number_of_fibers(i), ...
            'fiber_distribution', distribution, 'volume', volume, ...
            'position', position);
        
        fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, cells, ...
            fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
        
        % figure out brightest neurons
        cell_bright = max(fiber_mix, [], 1);
        [cell_bright_sort, idx] = sort(cell_bright, 'descend');
        
        % visbility of top neurons
        m = sort(fiber_mix(:, idx(1:top_neurons)), 'descend');
        if norm
            m = bsxfun(@rdivide, m, m(1, :));
        end
        
        result(i, :, :, j) = reshape(m(1:top_fibers, :), 1, top_fibers, top_neurons, 1);
    end
end

rng(old_rng);

%%

h = figure;
h.Position = h.Position .* [1 1 2 2];

for i = 1:length(number_of_fibers)
    subplot(2, 2, i);
    boxplot(reshape(squeeze(result(i, :, :, :)), top_fibers, [])' .* 100);
    xlabel('Fibers in descending order of fluence');
    if norm
        ylabel('% of highest fluence fiber');
    else
        ylabel('% of max');
    end
    title(sprintf('%d fibers', number_of_fibers(i)));
end

%r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig4-fibers.eps'); close;

%% figure 5: sample signals

close all;

% constants
number_of_fibers = 250; % good results: 100
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5]; % good results: [25 0 0; 0 25 0; 0 0 5]

old_rng = rng; rng(1);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
    'duration', 200, 'number_of_outputs', number_of_fibers, 'output_noise', 0, ...
    'params_fibers', {'fiber_distribution', distribution, 'angle_distribution', 0, 'volume', volume, 'position', position}, ...
    'params_cells', {'volume', volume}, 'figures', 2);
rng(old_rng);

%r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig5-signals.eps'); close;
%close all;

%% figure 6: auc

close all;

old_rng = rng; rng(0);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, 'duration', 200, 'number_of_outputs', 100, 'auc_threshold', 0.2:0.1:0.8);
rng(old_rng);

h = figure(7);
h.Position = [1000 200 880 812];
r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig6-auc.eps'); close;
close all;

%% figure ?: cell density vs fiber density
%explore({'mode', 'profile-rt', 'duration', 200, 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, 'params_fibers', {'angle_distribution', 0, 'fiber_distribution', [100 0 0; 0 100 0; 0 0 0], 'distribution', 'uniform'}}, ...
%    'number_of_outputs', 50:50:500, ...
%    'param_b_name', param_b_values, 5, true);
