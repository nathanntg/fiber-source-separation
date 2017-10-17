%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 4);
set(0, 'DefaultAxesFontSize', 24);

%% profile
% excitation
camera_profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
camera_profile_exc = sp_3d_to_2d(camera_profile_exc); % symmetric, way faster

im = [camera_profile_exc.volume(end:-1:2, :); camera_profile_exc.volume];
y = [-1 * camera_profile_exc.r(end:-1:2) camera_profile_exc.r];
figure;
imagesc(camera_profile_exc.z, y, log10(im), [-4 0]);
%title('Excitation profile');
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
xlim([-75 675]); ylim([-375 375]); axis square;
xticks([0 300 600]); yticks([-300 0 300]);
colorbar('Ticks', [-4 -2 0], 'TickLabels', {'10^{-4}', '10^{-2}', '10^0'});
colormap('jet');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/profile.eps'); close;

%% simulate: camera
% excitation
camera_profile_exc = sp_model('sensitivity-profile/camera-exc.mat');
camera_profile_exc = sp_3d_to_2d(camera_profile_exc); % symmetric, way faster

% generate mixing
camera_width = 64; camera_height = 48;
camera_number = camera_width * camera_height;
old_rng = rng(0);
[camera_mix, ~, camera_cells] = generate_camera_mixing(camera_width, camera_height, camera_profile_exc, 'figures', false);
rng(old_rng);

% useful values
camera_signal = max(camera_mix, [], 2);
camera_background = sum(camera_mix, 2) - camera_signal;


%% simulate: fibers
% excitation
fiber_profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
fiber_profile_exc = sp_3d_to_2d(fiber_profile_exc); % symmetric, way faster

% generate mixing
fiber_number = camera_number;
old_rng = rng(0);
[fiber_mix, ~, fiber_cells] = generate_realistic_mixing(fiber_number, fiber_profile_exc, ...
    'figures', false, 'fiber_distribution', [250 0 0; 0 250 0; 0 0 30]);
rng(old_rng);

% useful values
fiber_signal = max(fiber_mix, [], 2);
fiber_background = sum(fiber_mix, 2) - fiber_signal;

%% camera simulated image
% clearest neurons
[value, clearest_cell] = max(camera_mix, [], 2);
[~, ~, clearest_num] = unique(clearest_cell);
im = label2rgb(reshape(clearest_num, camera_height, camera_width), 'lines');

% figure w/o alpha
figure;
image(im);
axis xy; xticks([]); yticks([]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/camera-cells.eps'); close;

% figure w/ alpha
im_alpha = mat2gray(reshape(value, camera_height, camera_width));
im_composite = mat2gray(im) .* repmat(im_alpha, 1, 1, 3);
figure;
image(im_composite);
axis xy; xticks([]); yticks([]);

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/camera-cells-composite.eps'); close;

%% camera cbr
% make histograms
bins = 41;
mx = max(max(camera_signal ./ camera_background), max(fiber_signal ./ fiber_background));
edges = linspace(0, mx, bins);
camera_cbr = histcounts(camera_signal ./ camera_background, edges);
fiber_cbr = histcounts(fiber_signal ./ fiber_background, edges);

figure;
plot(edges(2:end), camera_cbr, edges(2:end), fiber_cbr);
xticks([0 0.05 0.1]); xlabel('Contrast to background ratio');
yticks([]); ylabel('Proportion'); % of pixels & fibers');
legend('Camera pixel', 'Fiber', 'Location', 'NorthEast');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/cbr.eps'); close;

%% depths
threshold = 0.05; % 0.05;

% camera
camera_depth = [];
camera_cbr = [];
camera_vis = [];
%cbr = camera_signal ./ camera_background;
%[~, clearest_cell] = max(camera_mix, [], 2);
%for j = unique(clearest_cell(:)')
for j = find(max(camera_mix, [], 1) > threshold)
    camera_depth = [camera_depth camera_cells(3, j)];
    camera_vis = [camera_vis max(camera_mix(:, j))];
    %camera_cbr = [camera_cbr max(cbr(clearest_cell == j))];
end

% fiber
fiber_depth = [];
fiber_cbr = [];
fiber_vis = [];
%cbr = fiber_signal ./ fiber_background;
%[~, clearest_cell] = max(fiber_mix, [], 2);
%for j = unique(clearest_cell(:)')
for j = find(max(fiber_mix, [], 1) > threshold)
    fiber_depth = [fiber_depth fiber_cells(3, j)];
    fiber_vis = [fiber_vis max(fiber_mix(:, j))];
    %fiber_cbr = [fiber_cbr max(cbr(clearest_cell == j))];
end

bins = 9;

% plot separately: camera
mn = min(camera_depth); mx = max(camera_depth);
[count, edges] = histcounts(camera_depth, bins);

figure;
h1 = plot(camera_depth, camera_vis * 100, '*');
xlabel('Depth [{\mu}]'); xlim([0 80]); xticks([0 40 80]);
ylabel('Visibility [%]'); yticks([0 50 100]);

yyaxis right;
h2 = plot(edges(2:end), count);
ylim([0 40]); yticks([0 20 40]);
set(gca, 'YColor', [0 0 0]);

legend([h1; h2], 'Neurons', 'Depth distribution', 'Location', 'NorthEast');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/depth-camera.eps'); close;

% plot separately: fiber
mn = min(fiber_depth); mx = max(fiber_depth);
[count, edges] = histcounts(fiber_depth, bins);

figure;
h1 = plot(fiber_depth, fiber_vis * 100, '*');
xlabel('Depth [{\mu}]'); xlim([180 300]); xticks([180 240 300]);
ylabel('Visibility [%]'); yticks([0 50 100]);
ylim([0 120]);

c = 200;
r = ylim();
line(c * [1 1], r, 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
text(c, r(2) - 2, 'Implant', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Rotation', 90);

yyaxis right;
h2 = plot(edges(2:end), count);
ylim([0 80]); yticks([0 40 80]);
set(gca, 'YColor', [0 0 0]);

legend([h1; h2], 'Neurons', 'Depth distribution', 'Location', 'NorthEast');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/depth-fiber.eps'); close;


