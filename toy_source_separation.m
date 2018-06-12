%% load profiles
% excitation
fiber_profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
fiber_profile_exc = sp_3d_to_2d(fiber_profile_exc); % symmetric, way faster

% emission
fiber_profile_emi = sp_model('sensitivity-profile/fiber-fluor.mat');
fiber_profile_emi = sp_3d_to_2d(fiber_profile_emi); % symmetric, way faster

%% setup constants
waveform = 'gcamp6f';
duration = 50; % seconds
mps = 100; % model time steps per second
sps = 20; % samples per second
spike_frequency = 0.4; % Hz (spikes per seconds)

volume = [100, 10, 100];
cells =  [60 50 40 10 20 30 40 50 60 70 80 90 ; ... x
           5  5  5  5  5  5  5  5  5  5  5  5 ; ... y
          50 50 50 80 80 80 80 80 80 80 80 80]; ... z
% cells =  [20 40 100 ; ... x
%            5  5  10 ; ... y
%           40 40 100]; ... z
fibers = [65 55 45; ... x
           5  5  5; ... y
          35 35 35]; ... z

%noise 
input_noise = 0; % either single value (0) or callback that generates noise
input_noise_type = 'scale';
output_noise = 0; % either single value (0) or callback that generates noise
output_noise_type = 'scale';

%% generate realistic mixing matrix
m = generate_realistic_rt_mixing(fibers, zeros(2, size(fibers, 2)), cells, fiber_profile_exc, fiber_profile_emi, ...
    'figures', false, 'stats', false);

figure; imagesc(m); title('Mixing matrix');

%% generate signals
% reset number generator
old_rng = rng; rng(1);

waveform_v = get_waveform(waveform, mps); % wave form in model samples
s = generate_inputs(size(cells, 2), spike_frequency, mps, sps, duration, waveform_v, [0.1 0.01], [1 2]);
%s(3, :) = 0;
duration_smp = size(s, 2); % duration in readout samples
s_noisy = add_noise(s, input_noise, input_noise_type);

% restore old number generator
rng(old_rng);

%% mix signals
x = m * s_noisy;
x_noisy = add_noise(x, output_noise, output_noise_type);

%% unmix
%[m_hat, s_hat] = unmix_als(x_noisy);
%[m_hat, s_hat] = unmix_nnmf_sparse(x_noisy, get_waveform(waveform, sps));
%[m_hat, s_hat] = unmix(x_noisy, get_waveform(waveform, sps));
%[s_hat, m_hat, w_hat] = fastica(x_noisy, 'g', 'skew', 'numOfIC', size(fibers, 2));
[s_hat, m_hat, w_hat] = unmix_nnica(x_noisy, get_waveform(waveform, sps));

%% visualize
h = figure; 
h.Position = [1 1 1 3] .* h.Position;

subplot(3, 1, 1);
plot_many(s');

subplot(3, 1, 2);
plot_many(x_noisy');

subplot(3, 1, 3);
plot_many(s_hat');

%% data to draw contours
c = load('sensitivity-profile/detail-exc.mat');
x = linspace(...
    double(c.config.voxel_min(1) - 1) * c.config.voxel_size(1),...
    double(c.config.voxel_max(1) - 1) * c.config.voxel_size(1),...
    double(c.config.voxel_count(1))) - c.config.initial_position(1);
z = linspace(...
    double(c.config.voxel_min(3) - 1) * c.config.voxel_size(3),...
    double(c.config.voxel_max(3) - 1) * c.config.voxel_size(3),...
    double(c.config.voxel_count(3))) - c.config.initial_position(3);

%% generate figure

n = 3;
coi = lines(n);
coo = [0.4 0.4 0.4];

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
line([18 38], [25 25], 'Color', 'black');
text(28, 23, '20 µm', 'FontSize', 22, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'black');

% panel b: example neural traces
subplot(2, 2, 2);
plot_many((1:size(s_noisy, 2)) ./ sps, s_noisy(1:n, :)');
xlim([0 20]);
title('Neural traces');

% panel c: recorded fiber signals
clrs = lines(n + size(fibers, 2));
subplot(2, 2, 3);
plot_many((1:size(x_noisy, 2)) ./ sps, x_noisy', clrs((n+1):end, :));
xlim([0 20]);
title('Signals via fibers');

% panel d: separated signals
rho = corr(s_noisy', s_hat');
[~, idx] = max(rho, [], 2);

subplot(2, 2, 4);
plot_many((1:size(s_hat, 2)) ./ sps, s_hat(idx(1:n), :)');
h = findobj(gca, 'Type', 'line');
% for i = 1:length(h)
%     h(i).LineStyle = ':';
% end
xlim([0 20]);
title('Separated signals');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Local/fig7-signals.eps'); close;
