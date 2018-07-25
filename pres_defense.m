%% setup
% sizing
set(groot, 'DefaultAxesLineWidth', 2);
set(groot, 'DefaultColorbarLineWidth', 2);
set(groot, 'DefaultLineLineWidth', 4);
set(groot, 'DefaultAxesFontSize', 24);
set(groot, 'DefaultAxesFontName', 'Helvetica');

% black background
set(groot, 'DefaultFigureColor', [0 0 0]);
set(groot, 'DefaultAxesXColor', [1 1 1] - 0.15);
set(groot, 'DefaultAxesYColor', [1 1 1] - 0.15);

% preserve background
set(groot, 'DefaultFigureInvertHardcopy', 'off');

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

print(gcf, '~/Local/profile.png', '-dpng', '-r300'); close;

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
        h = text(.98 * l, .98 * -l, sprintf('Neurons: %d', length(cell_bright_sort)), 'FontSize', 22, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', [1 1 1]);
        
        set(gca, 'Color', [0.2 0.2 0.2]);
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
    'Units', 'normalized', 'Color', [1 1 1]);
text(ax, spp2(1) + spp2(3) / 2, spp2(2) + spp2(4) + 0.02, sprintf('Fibers: %d', number_of_fibers(2)), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized', 'Color', [1 1 1]);
text(ax, spp1(1) - 0.07, spp1(2) + spp1(4) / 2, sprintf('Density: %dk per mm^3', cell_density(1) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized', 'Rotation', 90, 'Color', [1 1 1]);
text(ax, spp3(1) - 0.07, spp3(2) + spp3(4) / 2, sprintf('Density: %dk per mm^3', cell_density(2) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Units', 'normalized', 'Rotation', 90, 'Color', [1 1 1]);

rng(old_rng);

print(gcf, '~/Local/splay.png', '-dpng', '-r300'); close;

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
    title('Neurons per fiber', 'Color', [1 1 1]);
    
    set(gca, 'Color', [0.2 0.2 0.2]);
end

% titles
ax = axes('Position', [0 0 1 1], 'Visible', 'off');
for cd = 1:length(cell_density)
    spp = pos{cd};
    text(ax, spp(1) - 0.07, spp(2) + spp(4) / 2, sprintf('Density: %dk per mm^3', cell_density(cd) * 1000^3 / 1000), 'FontWeight', 'bold', 'FontSize', 24, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Units', 'normalized', 'Rotation', 90, 'Color', [1 1 1]);
end

rng(old_rng);

print(gcf, '~/Local/cells.png', '-dpng', '-r300'); close;

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

%% toy model

close all;

toy_source_separation;
close all;

%% recreate

n = 3;
coi = lines(n);
coo = [0.6 0.6 0.6];

h = figure('Renderer', 'painters');
h.Position = [1 1 2.1 2.1] .* h.Position;
    
% panel a: position of cells and fibers
subplot(2, 2, 1);
hold on;
scatter(cells(3, 1:n), cells(1, 1:n), 70, coi, 'filled'); % target cells
scatter(cells(3, (n+1):end), cells(1, (n+1):end), 50, coo, 'filled'); % background cells
for i = 1:size(fibers, 2)
    line([0 fibers(3, i)], [1 1] .* fibers(1, i), 'Color', coo, 'LineWidth', 2);
    contour(fibers(3, i) + z .* 1000,fibers(1, i) + x .* 1000, log10(max(mat2gray(squeeze(c.I(:, 25, :))), 1e-5)), [-3 -2 -1], ':', 'LineWidth', 1.5);
    colormap('jet'); caxis([-4 0]);
end

hold off;
axis square;
xlim([15 85]); xticks([]); xlabel('z [µm]');
ylim([15 85]); yticks([]); ylabel('x [µm]');

% scale bar
line([18 38], [25 25], 'Color', 'white');
text(28, 23, '20 µm', 'FontSize', 22, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'white');
set(gca, 'Color', [0.2 0.2 0.2]);

% panel b: example neural traces
subplot(2, 2, 2);
plot_many((1:size(s_noisy, 2)) ./ sps, s_noisy(1:n, :)');
xlim([0 20]); xticks([0 10 20]); xlabel('Time [s]');
title('Neural traces', 'Color', [1 1 1]);
set(gca, 'Color', [0.2 0.2 0.2]);

% panel c: recorded fiber signals
clrs = lines(n + size(fibers, 2));
subplot(2, 2, 3);
plot_many((1:size(x_noisy, 2)) ./ sps, x_noisy', clrs((n+1):end, :));
xlim([0 20]); xticks([0 10 20]); xlabel('Time [s]');
title('Signals via fibers', 'Color', [1 1 1]);
set(gca, 'Color', [0.2 0.2 0.2]);

% panel d: separated signals
rho = corr(s_noisy', s_hat');
[~, idx] = max(rho, [], 2);

subplot(2, 2, 4);
plot_many((1:size(s_hat, 2)) ./ sps, s_hat(idx(1:n), :)');
h = findobj(gca, 'Type', 'line');
% for i = 1:length(h)
%     h(i).LineStyle = ':';
% end
xlim([0 20]); xticks([0 10 20]); xlabel('Time [s]');
title('Separated signals', 'Color', [1 1 1]);
set(gca, 'Color', [0.2 0.2 0.2]);

print(gcf, '~/Local/separate.png', '-dpng', '-r300'); close;

%% compare sensitivity profile

h = figure();
h.Position = h.Position .* [1 1 2 1.1];

% draw image
subplot(1, 2, 1);
photo = imread('/Users/nathan/Documents/School/BU/Gardner Lab/Fiber/Modeling Paper/profile/fitc.jpg');
photo = permute(photo, [2 1 3]);
image(photo);
axis xy; xlabel('z [µm]'); ylabel('x [µm]');
xticks([]); yticks([]); axis square;

% 70px/6um
xlim(880 + [-500 500]); ylim(1370 + [-500 500]);

% scale bar
line([420 653], [980 980], 'Color', 'white');
text(536.5, 960, '20 µm', 'FontSize', 22, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'white');

% profile in water
detail_profile_exc_water = sp_model('sensitivity-profile/detail-exc-water.mat');
detail_profile_exc_water = sp_3d_to_2d(detail_profile_exc_water);

ax2 = subplot(1, 2, 2);
im = [detail_profile_exc_water.volume(end:-1:2, :); detail_profile_exc_water.volume];
y = [-1 * detail_profile_exc_water.r(end:-1:2) detail_profile_exc_water.r];
imagesc(detail_profile_exc_water.z, y, im, [0 1]);
%title('Excitation profile');
axis xy; xlabel('z [µm]'); ylabel('x [µm]');
xlim([-20 65]); ylim([-42.5 42.5]); axis square;
xticks([]); yticks([]);
colormap('jet'); caxis([0 1]);

% colorbar position
spp = get(ax2, 'Position');
colorbar('Ticks', [0 0.5 1], ...
    'Position', [spp(1) + spp(3) + 0.02, spp(2), 0.01, spp(4)]);

print(gcf, '~/Local/compare-sensitivity.png', '-dpng', '-r300'); close;

%r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig2-profile.eps'); close;

%% Robustness

param_fiber_count = [100 250]; % 500];
param_spike_frequency = [0.2 0.4 0.8];

load('state.mat');

% actual plotting
figure;
h = bar(results_mean(:, 2:3));

for i = 1:length(param_fiber_count)
    fiber_count = param_fiber_count(i);
    
    idx_start = (i - 1) * length(param_spike_frequency) + 1;
    idx_end = (i - 1) * length(param_spike_frequency) + length(param_spike_frequency);
    
    line([idx_start - 0.4 idx_end + 0.4], [1.11 1.11], 'Color', 'white', 'LineWidth', 1);
    text((idx_start + idx_end) / 2, 1.11, sprintf('%d fibers', fiber_count), 'Color', 'white', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    for j = 1:length(param_spike_frequency)
        spike_frequency = param_spike_frequency(j);
        
        idx = (i - 1) * length(param_spike_frequency) + j;
        
        line([idx - 0.4 idx + 0.4], [1.01 1.01], 'Color', 'white', 'LineWidth', 1);
        text(idx, 1.01, sprintf('%.1f Hz', spike_frequency), 'Color', 'white', 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        line([idx - 0.14 idx - 0.14], [results_mean(idx, 2) - results_std(idx, 2) / 2 results_mean(idx, 2) + results_std(idx, 2) / 2], 'Color', 'white', 'LineWidth', 1);
        %line([idx idx], [results_mean(idx, 2) - results_std(idx, 2) / 2 results_mean(idx, 2) + results_std(idx, 2) / 2], 'Color', 'white', 'LineWidth', 1);
        line([idx + 0.14 idx + 0.14], [results_mean(idx, 3) - results_std(idx, 3) / 2 results_mean(idx, 3) + results_std(idx, 3) / 2], 'Color', 'white', 'LineWidth', 1);
    end
end

xticks([]);
ylabel('Accurately matched [%]'); ylim([0 1.2]); yticks([0 0.5 1.0]);
set(gca, 'Color', [0.2 0.2 0.2]);

h2 = legend('Raw', 'NN-ICA', 'Location', 'SouthOutside', 'Orientation', 'horizontal');
set(h2, 'Color', [0.2 0.2 0.2]);
set(h2, 'TextColor', [1 1 1]);

print(gcf, '~/Local/robust.png', '-dpng', '-r150'); close;

%% ROC curve (Backup)

% idea: https://stats.stackexchange.com/questions/186337/average-roc-for-repeated-10-fold-cross-validation-with-probability-estimates

iters = 5; % number of iterations
thresholds = [0.2 0.25 0.3 0.5 0.7 0.9];
colors = lines(length(thresholds));
base_fpr = linspace(0, 1, 101);
tprs = zeros(length(thresholds), length(base_fpr), iters);
nums = zeros(length(thresholds), iters);
aucs = zeros(length(thresholds), iters);

% create figure
h = figure;
h.Position = [1 1 880 812];

hold on;

old_rng = rng; rng(0);
for i = 1:iters
    [~, ~, auc] = simulate_source_separation('mode', 'profile-rt', ...
        'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
        'duration', 200, 'number_of_outputs', 100, 'auc_threshold', thresholds, ...
        'output_noise', 0, 'g', 'nnica', 'params_cells', {'cell_density', 0.00025}, ...
        'sss_mode', 'auc', 'figures', false);
    
    for j = 1:length(thresholds)
        fpr = auc(j).fpr;
        tpr = auc(j).tpr;
        
        % draw individual trials
        %plot(fpr, tpr, 'LineWidth', 1, 'Color', colors(j, :));
        
        % other values
        nums(j, i) = auc(j).number;
        aucs(j, i) = auc(j).auc;
        
        % unique
        [fpr, idx] = unique(fpr, 'last');
        tpr = tpr(idx);
        if ~isempty(tpr)
            tprs(j, :, i) = interp1(fpr, tpr, base_fpr, 'linear', 'extrap');
            tprs(j, 1, i) = 0;
        else
            tprs(j, 2:end, i) = 1;
        end
    end
end
rng(old_rng);

% draw mean lines
l = cell(1, length(thresholds));
h = zeros(1, length(thresholds));
for j = 1:length(thresholds)
    % mean
    mn = mean(tprs(j, :, :), 3);
    st = std(tprs(j, :, :), 0, 3);
    
    % draw patch
    %patch_x = [base_fpr fliplr(base_fpr)];
    %patch_y = [max(mn - st, 0) fliplr(min(mn + st, 1))];
    %patch(patch_x, patch_y, 1, 'FaceColor', colors(j, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);

    % draw line
    h(j) = plot(base_fpr, mn, 'Color', colors(j, :));
    auc = trapz(base_fpr, mn);
    
    % make legend entry
    l{j} = sprintf('r^2 > %.2f; AUC: %.3f; N: %.1f', thresholds(j), mean(aucs(j, :)), mean(nums(j, :)));
end

hold off;

xlabel('False Positive Rate'); xlim([0 1]); xticks([0 0.5 1.0]);
ylabel('True Positive Rate'); ylim([0 1]); yticks([0 0.5 1.0]);
set(gca, 'Color', [0.2 0.2 0.2]);

h2 = legend(h, l, 'Location', 'SouthEast');
set(h2, 'Color', [0.2 0.2 0.2]);
set(h2, 'TextColor', [1 1 1]);

%print(gcf, '~/Local/roc.png', '-dpng', '-r150'); close;
