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

%% figure 2: sensitivity profile
h = figure();
h.Position = h.Position .* [1 1 3 1.1];

% draw image
subplot(1, 3, 1);
photo = imread('/Users/nathan/Documents/School/BU/Gardner Lab/Fiber/Modeling Paper/profile/fitc.jpg');
photo = permute(photo, [2 1 3]);
photo = rgb2gray(im2double(photo));
imagesc(photo);
axis xy; xlabel('z [µm]'); ylabel('x [µm]');
xticks([]); yticks([]); axis square;

% 70px/6um
xlim(880 + [-500 500]); ylim(1370 + [-500 500]);

% scale bar
line([420 653], [980 980], 'Color', 'white');
text(536.5, 960, '20 µm', 'FontSize', 22, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'white');

% profile in water
detail_profile_exc_water = sp_model('sensitivity-profile/detail-exc-water.mat');
%detail_profile_exc_water = sp_3d_to_2d(detail_profile_exc_water);

subplot(1, 3, 2);
im = squeeze(sum(detail_profile_exc_water.volume, 2));
im = im ./ max(im(:));
y = detail_profile_exc_water.x;
imagesc(detail_profile_exc_water.z, y, im, [0 1]);
%title('Excitation profile');
axis xy; xlabel('z [µm]'); ylabel('x [µm]');
xlim([-20 65]); ylim([-42.5 42.5]); axis square;
xticks([]); yticks([]);
colormap('jet'); caxis([0 1]);

% profile in water
detail_profile_exc = sp_model('sensitivity-profile/detail-exc.mat');
%detail_profile_exc = sp_3d_to_2d(detail_profile_exc);

ax3 = subplot(1, 3, 3);
im = squeeze(sum(detail_profile_exc.volume, 2));
im = im ./ max(im(:));
y = detail_profile_exc.x;
imagesc(detail_profile_exc.z, y, im, [0 1]);
%title('Excitation profile');
axis xy; xlabel('z [µm]'); ylabel('x [µm]');
xlim([-20 65]); ylim([-42.5 42.5]); axis square;
xticks([]); yticks([]);
colormap('jet'); caxis([0 1]);

% colorbar position
spp = get(ax3, 'Position');
colorbar('Ticks', [0 0.5 1], ...
    'Position', [spp(1) + spp(3) + 0.02, spp(2), 0.01, spp(4)]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig2-profile.eps'); close;

%% figure 3: fiber distribution

% constants
number_of_fibers = [250 500];
cell_density = [0.00025 0.0005];
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

% open window
h = figure('Renderer', 'painters');
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
        xlim([-l l]); xticks([-500 0 500]); xlabel('x [µm]');
        ylim([-l l]); yticks([-500 0 500]); ylabel('y [µm]');
        
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

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig3-distribution.eps'); close;

%% figure 4: number of cells seen

% constants
iters = 5;
thresholds = [0.01 0.02 0.05 0.1];
number_of_fibers = [50 100 150 200 500 750 1000 1500 2000];
cell_density = [0.00025 0.0005];
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

% plot it
h = figure('Renderer', 'painters');
h.Position(3) = 2 * h.Position(3);
h.Position(4) = length(cell_density) * h.Position(4);

old_rng = rng; rng(0);

pos = {};

for cd = 1:length(cell_density)
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
    ylim([10 99000]);
    xlim(number_of_fibers([1 end]));

    l = cell(1, length(thresholds));
    for i = 1:length(thresholds)
        l{i} = sprintf('\\geq %.1f%%', thresholds(i)*100);
    end
    legend(h, l{:}, 'Location', 'NorthWest');

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
    
    ylim([0.2 500]); yticks([1 10 100]);

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

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig4-cells.eps'); close;

%% figure 5: neuron strengths per fiber

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
cells = generate_cells('volume', volume, 'cell_density', 0.00025);

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

h = figure('Renderer', 'painters');
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

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig5-neurons.eps'); close;

%% figure 6: number of fibers per cell

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
cells = generate_cells('volume', volume, 'cell_density', 0.00025);

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

h = figure('Renderer', 'painters');
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

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig6-fibers.eps'); close;

%% figure 7: sample signals

close all;
toy_source_separation;

%% figure 8: compare performance

% two parameters to explore: firing rate, number of fibers

iters = 3; % number of iterations per parameter set
param_fiber_count = [100 250]; % 500];
param_spike_frequency = [0.2 0.4 0.8];

% iters = 1;
% param_fiber_count = [100 250];
% param_spike_frequency = 0.4;

threshold = 0.6;

% used cache results, since relatively slow
if false || ~exist('state.mat', 'file')
    results_mean = zeros(length(param_spike_frequency) * length(param_fiber_count), 3);
    results_std = zeros(length(param_spike_frequency) * length(param_fiber_count), 3);

    old_rng = rng; rng(0);
    for i = 1:length(param_fiber_count)
        for j = 1:length(param_spike_frequency)
            fiber_count = param_fiber_count(i);
            spike_frequency = param_spike_frequency(j);

            fprintf('Simulating %d fibers with spike frequency of %.1f\n', fiber_count, spike_frequency);

            values = zeros(iters, 3);

            for k = 1:iters
                % control
                scores = simulate_source_separation('mode', 'profile-rt', ...
                    'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
                    'duration', 200, 'number_of_outputs', fiber_count, 'spike_frequency', spike_frequency, ...
                    'output_noise', 0, 'g', 'control', 'params_cells', {'cell_density', 0.00025}, ...
                    'figures', false);

                values(k, 1) = sum(scores > threshold) ./ fiber_count;

                % no source separation
                scores = simulate_source_separation('mode', 'profile-rt', ...
                    'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
                    'duration', 200, 'number_of_outputs', fiber_count, 'spike_frequency', spike_frequency, ...
                    'output_noise', 0, 'g', 'none', 'params_cells', {'cell_density', 0.00025}, ...
                    'figures', false);

                values(k, 2) = sum(scores > threshold) ./ fiber_count;

                % regular source separation
                scores = simulate_source_separation('mode', 'profile-rt', ...
                    'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
                    'duration', 200, 'number_of_outputs', fiber_count, 'spike_frequency', spike_frequency, ...
                    'output_noise', 0, 'g', 'nnica', 'params_cells', {'cell_density', 0.00025}, ...
                    'figures', false);

                values(k, 3) = sum(scores > threshold) ./ fiber_count;
            end

            % incremental logging of sorts
            disp([i j]);
            disp(values);

            % average
            mn = mean(values, 1);
            st = std(values, 0, 1);

            % store
            idx = (i - 1) * length(param_spike_frequency) + j;
            results_mean(idx, :) = mn;
            results_std(idx, :) = st;
        end
    end
    rng(old_rng);

    % save the data, it's sloooooooooow
    save('state.mat', '-v7.3', 'results_mean', 'results_std');
else
    load('state.mat');
end

% actual plotting
figure;
h = bar(results_mean);

for i = 1:length(param_fiber_count)
    fiber_count = param_fiber_count(i);
    
    idx_start = (i - 1) * length(param_spike_frequency) + 1;
    idx_end = (i - 1) * length(param_spike_frequency) + length(param_spike_frequency);
    
    line([idx_start - 0.4 idx_end + 0.4], [1.11 1.11], 'Color', 'black', 'LineWidth', 1);
    text((idx_start + idx_end) / 2, 1.11, sprintf('%d fibers', fiber_count), 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    for j = 1:length(param_spike_frequency)
        spike_frequency = param_spike_frequency(j);
        
        idx = (i - 1) * length(param_spike_frequency) + j;
        
        line([idx - 0.4 idx + 0.4], [1.01 1.01], 'Color', 'black', 'LineWidth', 1);
        text(idx, 1.01, sprintf('%.1f Hz', spike_frequency), 'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        line([idx - 0.22 idx - 0.22], [results_mean(idx, 1) - results_std(idx, 1) / 2 results_mean(idx, 1) + results_std(idx, 1) / 2], 'Color', 'black', 'LineWidth', 1);
        line([idx idx], [results_mean(idx, 2) - results_std(idx, 2) / 2 results_mean(idx, 2) + results_std(idx, 2) / 2], 'Color', 'black', 'LineWidth', 1);
        line([idx + 0.22 idx + 0.22], [results_mean(idx, 3) - results_std(idx, 3) / 2 results_mean(idx, 3) + results_std(idx, 3) / 2], 'Color', 'black', 'LineWidth', 1);
    end
end

xticks([]);
ylabel('Accurately matched [%]'); ylim([0 1.2]); yticks([0 0.5 1.0]);
legend('Control (random)', 'Control (raw)', 'NN-ICA', 'Location', 'SouthOutside', 'Orientation', 'horizontal');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig8-robust.eps'); close;

%% figure 9: auc
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

legend(h, l, 'Location', 'SouthEast');
xlabel('False Positive Rate'); xlim([0 1]); xticks([0 0.5 1.0]);
ylabel('True Positive Rate'); ylim([0 1]); yticks([0 0.5 1.0]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig9-auc.eps'); close;

%% figure ?: auc control 1: compare with random traces

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
        'output_noise', 0, 'g', 'unmix', 'params_cells', {'cell_density', 0.00025}, ...
        'sss_mode', 'auc_control', 'figures', false);
    
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
        tprs(j, :, i) = interp1(fpr, tpr, base_fpr, 'linear', 'extrap');
        tprs(j, 1, i) = 0;
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

legend(h, l, 'Location', 'SouthEast');
xlabel('False Positive Rate'); xlim([0 1]); xticks([0 0.5 1.0]);
ylabel('True Positive Rate'); ylim([0 1]); yticks([0 0.5 1.0]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig-auc-control1.eps'); close;

%% figure ?: auc control 2; no unmixing

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
        'output_noise', 0, 'g', 'none', 'params_cells', {'cell_density', 0.00025}, ...
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
        tprs(j, :, i) = interp1(fpr, tpr, base_fpr, 'linear', 'extrap');
        tprs(j, 1, i) = 0;
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

legend(h, l, 'Location', 'SouthEast');
xlabel('False Positive Rate'); xlim([0 1]); xticks([0 0.5 1.0]);
ylabel('True Positive Rate'); ylim([0 1]); yticks([0 0.5 1.0]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig-auc-control2.eps'); close;

%% figure ?: cell density vs fiber density
%explore({'mode', 'profile-rt', 'duration', 200, 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, 'params_fibers', {'angle_distribution', 0, 'fiber_distribution', [100 0 0; 0 100 0; 0 0 0], 'distribution', 'uniform'}}, ...
%    'number_of_outputs', 50:50:500, ...
%    'param_b_name', param_b_values, 5, true);


%% figure ?: multiplexing feasibility

opt_fiber = 1e-3;
opt_neuron = 1e-4;
number_of_fibers = 1000;
volume = [1200; 1200; 400];
position = [0.5; 0.5; 0.25];
distribution = [150 0 0; 0 150 0; 0 0 5];

for i = 1:length(number_of_fibers)
    % generate fibers
    [fibers, fiber_angles] = generate_fibers(number_of_fibers(i));
    
    % generate cells
    cells = generate_cells('cell_density', 0.00025);
    
    % normal
    %m_alt = generate_realistic_rt_mixing(fibers, fiber_angles, cells, fiber_profile_exc, fiber_profile_emi, 'figures', false, 'stats', false);
    %m_alt = m_alt(:, any(m_alt > opt_neuron, 1));
    
    % mixing matrix
    m = generate_multiplex_mixing(fibers, fiber_angles, cells, fiber_profile_exc, fiber_profile_emi, number_of_fibers(i) / 2, false, 'opt_fiber', opt_fiber, 'opt_neuron', opt_neuron, 'as_single', true);
end

disp(size(m));

figure;
semilogy(svd(m * m'));
title('1000 fibers; 500 step multiplexing');


