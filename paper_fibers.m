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

%% generate profile
im = [fiber_profile_exc.volume(end:-1:2, :); fiber_profile_exc.volume];
y = [-1 * fiber_profile_exc.r(end:-1:2) fiber_profile_exc.r];
figure;
imagesc(fiber_profile_exc.z, y, log10(im), [-4 0]);
%title('Excitation profile');
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
xlim([-75 675]); ylim([-375 375]); axis square;
xticks([0 300 600]); yticks([-300 0 300]);
colorbar('Ticks', [-4 -2 0], 'TickLabels', {'10^{-4}', '10^{-2}', '10^0'});
colormap('jet');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/profile.eps'); close;

%% simulate: fibers
% generate mixing
fiber_number = 500;
old_rng = rng(0);
volume = [1200; 1200; 1200];
[fibers, fiber_angles] = generate_fibers(fiber_number, ...
    'position', [0.5; 0.5; 100], 'volume', volume);
fiber_cells = generate_cells('volume', volume);
fiber_mix = generate_realistic_rt_mixing(fibers, fiber_angles, fiber_cells, ...
    fiber_profile_exc, fiber_profile_emi, ...
    'normalize', true, 'figures', false);
rng(old_rng);

% useful values
fiber_signal = max(fiber_mix, [], 2);
fiber_background = sum(fiber_mix, 2) - fiber_signal;

%% simulate: chr2
old_rng = rng(0);
[~, chr2_exc_cells, chr2_exc_fibers] = generate_realistic_rt_mixing(fibers, fiber_angles, fiber_cells, ...
    fiber_profile_exc, fiber_profile_emi, 'normalize', true, 'figures', false, 'stats', false);
chr2_exc_cells = chr2_exc_cells ./ max(chr2_exc_fibers); % normalize
rng(old_rng);

%% depths
threshold = 0.01; % 0.05;

% fiber
fiber_depth = [];
fiber_vis = [];
%cbr = fiber_signal ./ fiber_background;
%[~, clearest_cell] = max(fiber_mix, [], 2);
%for j = unique(clearest_cell(:)')
for j = find(max(fiber_mix, [], 1) > threshold)
    fiber_depth = [fiber_depth fiber_cells(3, j)];
    fiber_vis = [fiber_vis max(fiber_mix(:, j))];
    %fiber_cbr = [fiber_cbr max(cbr(clearest_cell == j))];
end

% chr2
chr2_depth = [];
chr2_vis = [];
for j = find(chr2_exc_cells > threshold)
    chr2_depth = [chr2_depth fiber_cells(3, j)];
    chr2_vis = [chr2_vis chr2_exc_cells(j)];
end


h = figure;
h.Position = [h.Position(1) h.Position(2) h.Position(3) * 2.1 h.Position(4) * 1.1];

% plot fiber
bins = 9;
hs = subplot(1, 2, 2);

[count, edges] = histcounts(fiber_depth, bins);

c = 100;
h1 = plot(fiber_depth - c, fiber_vis * 100, '.', 'MarkerSize', 25);
xlabel('Depth [{\mu}]'); xlim([-35 125]); xticks([0 40 80 120]);
ylabel('Fluorescence yield [% of max]'); yticks([0 50 100]);
ylim([0 100]);

r = ylim();
line([0 0], r, 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
text(0, r(2) - 2, 'Implant', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Rotation', 90);

yyaxis right;
h2 = plot(edges(2:end) - c, count);
ylim([0 600]); yticks([0 250 500]);
set(gca, 'YColor', [0 0 0]);
h = ylabel('Number of neurons');
h.Rotation = -90;
h.VerticalAlignment = 'bottom';

legend([h1; h2], 'Neurons', 'Depth distribution', 'Location', 'NorthEast');

hs.Position(1) = hs.Position(1) + 0.03;
hs.Position(2) = hs.Position(2) + 0.04;

% plot chr2
bins = 19;
hs = subplot(1, 2, 1);

[count, edges] = histcounts(chr2_depth, bins);

c = 100;
h1 = plot(chr2_depth - c, chr2_vis * 100, '.', 'MarkerSize', 25);
xlabel('Depth [{\mu}]'); xlim([-50 800]); xticks([0 200 400 600 800]);
ylabel('Excitation [% of max]'); yticks([0 50 100]);
ylim([0 100]);

r = ylim();
line([0 0], r, 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
text(0, r(2) - 2, 'Implant', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Rotation', 90);

yyaxis right;
h2 = plot(edges(2:end) - c, count);
ylim([0 12500]); yticks([0 5000 10000]);
set(gca, 'YColor', [0 0 0]);
h = ylabel('Number of neurons');
h.Rotation = -90;
h.VerticalAlignment = 'bottom';

legend([h1; h2], 'Neurons', 'Depth distribution', 'Location', 'NorthEast');

hs.Position(1) = hs.Position(1) - 0.05;
hs.Position(2) = hs.Position(2) + 0.04;

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/depth-fiber.eps'); close;
