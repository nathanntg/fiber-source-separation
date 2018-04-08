%% setup
% sizing
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultLineLineWidth', 4);
set(groot, 'DefaultAxesFontSize', 24);

%% load profiles
% excitation
fiber_profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
fiber_profile_exc = sp_3d_to_2d(fiber_profile_exc); % symmetric, way faster

% emission
fiber_profile_emi = sp_model('sensitivity-profile/fiber-fluor.mat');
fiber_profile_emi = sp_3d_to_2d(fiber_profile_emi); % symmetric, way faster

%% profile
im = [fiber_profile_exc.volume(end:-1:2, :); fiber_profile_exc.volume];
y = [-1 * fiber_profile_exc.r(end:-1:2) fiber_profile_exc.r];
figure;
imagesc(fiber_profile_exc.z, y, log10(im), [-4 0]);
%title('Excitation profile', 'Color', [1 1 1]);
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
xlim([-75 675]); ylim([-375 375]); axis square;
xticks([0 300 600]); yticks([-300 0 300]);
colorbar('Ticks', [-4 -2 0], 'TickLabels', {'10^{-4}', '10^{-2}', '10^0'});
colormap('jet');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/profile.eps'); close;

%% fiber distribution

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
        text(.98 * l, .98 * -l, sprintf('Neurons: %d', length(cell_bright_sort)), 'FontSize', 22, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', [0 0 0]);
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

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/splay.eps'); close;

%% number of cells seen

% constants
number_of_fibers = [50 100 150 200 500 750 1000 1500 2000];
cell_density = [0.000275 0.00078];
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

% plot it
h = figure;
h.Position(3) = 2 * h.Position(3);
h.Position(4) = length(cell_density) * h.Position(4);

old_rng = rng; rng(0);

pos = {};

for cd = 1:length(cell_density)
    iters = 5;
    thresholds = [0.01 0.02 0.05 0.1];

    result_single = zeros(length(number_of_fibers), length(thresholds), iters);
    result_multi = zeros(length(number_of_fibers), length(thresholds), iters);
    result_per_fiber = zeros(length(number_of_fibers), length(thresholds), iters);
    result_per_fiber_multi = zeros(length(number_of_fibers), length(thresholds), iters);

    % generate fibers
    cells = generate_cells('volume', volume, 'cell_density', cell_density(cd));

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
                cell_vis_to_multi = sum(fiber_mix >= thresholds(k), 1) >= 2;
                result_single(i, k, j) = sum(cell_bright >= thresholds(k));
                result_multi(i, k, j) = sum(cell_vis_to_multi);
                result_per_fiber(i, k, j) = mean(sum(fiber_mix >= thresholds(k), 2));
                result_per_fiber_multi(i, k, j) = mean(sum(fiber_mix(:, cell_vis_to_multi) >= thresholds(k), 2));
            end
        end
    end

    hs = subplot(length(cell_density), 2, 2 * (cd - 1) + 1);
    pos{end + 1} = get(hs, 'Position');

    mn = mean(result_single, 3);
    % st = std(result_single, 0, 3);
    set(gca, 'ColorOrderIndex', 1);
    h = loglog(number_of_fibers, mn);
    % hold on;
    % plot(number_of_fibers, mn + st, 'LineWidth', 1);
    % plot(number_of_fibers, mn - st, 'LineWidth', 1);
    % hold off;
    title('Neurons', 'Color', [1 1 1]);
    xlabel('Number of fibers');
    ylabel('Neurons above threshold');

    hold on;
    mn = mean(result_multi, 3);
    % st = std(result_single, 0, 3);
    set(gca, 'ColorOrderIndex', 1);
    loglog(number_of_fibers, mn, ':');
    hold off;

    % adjust dimensions
    ylim([10 99000]);
    xlim(number_of_fibers([1 end]));
    
    set(gca, 'Color', [0.2 0.2 0.2]);

    l = cell(1, length(thresholds));
    for i = 1:length(thresholds)
        l{i} = sprintf('\\geq %.1f%%', thresholds(i)*100);
    end
    h2 = legend(h, l{:}, 'Location', 'NorthWest');
    set(h2, 'Color', [0.2 0.2 0.2]);
    set(h2, 'TextColor', [1 1 1]);

    % neurons per fiber
    subplot(length(cell_density), 2, 2 * (cd - 1) + 2);

    mn = mean(result_per_fiber, 3);
    set(gca, 'ColorOrderIndex', 1);
    loglog(number_of_fibers, mn);

    hold on;
    mn = mean(result_per_fiber_multi, 3);
    % st = std(result_single, 0, 3);
    set(gca, 'ColorOrderIndex', 1);
    loglog(number_of_fibers, mn, ':');
    hold off;

    xlim(number_of_fibers([1 end]));
    
    ylim([0.2 500]);

    xlabel('Number of fibers');
    ylabel('Neurons above threshold');
    title('Neurons per fiber');
end

% titles
ax = axes('Position', [0 0 1 1], 'Visible', 'off');
for cd = 1:length(cell_density)
    spp = pos{cd};
    text(ax, spp(1) - 0.07, spp(2) + spp(4) / 2, sprintf('Density: %dk per mm^3', cell_density(cd) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Units', 'normalized', 'Rotation', 90);
end

rng(old_rng);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/cells.eps'); close;


%% neuron strengths per fiber (UNUSED)

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
    h2 = boxplot(reshape(permute(result{i}, [2 1 3]), top_neurons, [])' .* 100);
    set(h2, 'LineWidth', 2);
    h3 = findobj(gca, 'tag', 'Outliers');
    for iH = 1:length(h3)
        h3(iH).MarkerEdgeColor = [0    0.4470    0.7410];
    end
    xlabel('Neurons in descending order of fluence');
    if norm
        ylabel('% of brightest');
    else
        ylabel('% of max');
    end
    title(sprintf('%d fibers', number_of_fibers(i)), 'Color', [1 1 1]);
    set(gca, 'Color', [0.2 0.2 0.2]);
end

%r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig3-neurons.eps'); close;

%% number of fibers per cell

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

h = figure;
h.Position = h.Position .* [1 1 2 2];

for i = 1:length(number_of_fibers)
    subplot(2, 2, i);
    h2 = boxplot(reshape(squeeze(result(i, :, :, :)), top_fibers, [])' .* 100, 'Color', [0.8500    0.3250    0.0980]);
    set(h2, 'LineWidth', 2);
    h3 = findobj(gca, 'tag', 'Outliers');
    for iH = 1:length(h3)
        h3(iH).MarkerEdgeColor = [0    0.4470    0.7410];
    end
    xlabel('Fibers in descending order of fluence');
    if norm
        ylabel('% of highest fluence fiber');
    else
        ylabel('% of max');
    end
    title(sprintf('%d fibers', number_of_fibers(i)), 'Color', [1 1 1]);
    set(gca, 'Color', [0.2 0.2 0.2]);
end

print(gcf, '~/Local/fibers.png', '-dpng', '-r300'); close;

%% sample signals (UNTESTED)

close all;

% constants
number_of_fibers = 250; % good results: 100
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5]; % good results: [25 0 0; 0 25 0; 0 0 5]

old_rng = rng; rng(1);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
    'duration', 200, 'number_of_outputs', number_of_fibers, ...
    'params_fibers', {'fiber_distribution', distribution, 'angle_distribution', 0, 'volume', volume, 'position', position}, ...
    'params_cells', {'volume', volume}, 'figures', 2, 'g', 'unmix');
rng(old_rng);

h = subplot(3, 1, 1);
set(h, 'Color', [0.2 0.2 0.2]); % background color
h2 = findobj(gcf, 'tag', 'Legend');
set(h2, 'Color', [0.2 0.2 0.2]);
set(h2, 'TextColor', [1 1 1]);

h = subplot(3, 1, 2);
set(h, 'Color', [0.2 0.2 0.2]); % background color
h2 = findobj(gcf, 'tag', 'Legend');
set(h2, 'Color', [0.2 0.2 0.2]);
set(h2, 'TextColor', [1 1 1]);

h = subplot(3, 1, 3);
set(h, 'Color', [0.2 0.2 0.2]); % background color

%print(gcf, '~/Local/separate.png', '-dpng', '-r300'); close;
%close all;
