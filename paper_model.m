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
number_of_fibers = 100;
volume = [1000; 1000; 400];
position = [0.5; 0.5; 0.25];
distribution = [125 0 0; 0 125 0; 0 0 5];

% generate fibers
old_rng = rng; rng(1);
[fibers, fiber_angles] = generate_fibers(number_of_fibers, ...
    'fiber_distribution', distribution, 'volume', volume, ...
    'position', position);
cells = generate_cells('volume', volume);
rng(old_rng);

% calculate center
center = volume .* position;

% open window
figure;
plot(fibers(1, :) - center(1), fibers(2, :) - center(2), 'g.', 'MarkerSize', 18);
axis square;
l = max(max(abs(fibers(1, :) - center(1))), max(abs(fibers(2, :) - center(2))));
l = 400; % ceil(l / 100) * 100;
xlim([-l l]); xticks([-l 0 l]); xlabel('x [{\mu}]');
ylim([-l l]); yticks([-l 0 l]); ylabel('y [{\mu}]');

r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig1-distribution.eps'); close;

%% figure 2: sample signals

close all;

old_rng = rng; rng(0);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', fiber_profile_exc, 'profile_fluor', fiber_profile_emi, ...
    'duration', 200, 'number_of_outputs', 100, 'output_noise', @(n) normrnd(0, 0.05, 1, n), ...
    'params_fibers', {'fiber_distribution', [25 0 0; 0 25 0; 0 0 5]});
rng(old_rng);

h = figure(4);
h.Position = [933 105 1200 1140];
r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig2-recovered.eps'); close;

h = figure(2);
h.Position = [1000 200 900 1140];
r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig3-fiber.eps'); close;

close all;

%% figure 3: auc
old_rng = rng; rng(0);
simulate_source_separation('mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 200, 'number_of_outputs', 100, 'auc_threshold', 0.2:0.1:0.8);
rng(old_rng);

h = figure(7);
h.Position = [1000 200 880 812];
r = get(gcf, 'renderer'); print(gcf, '-depsc2', ['-' r], '~/Desktop/fig4-auc.eps'); close;
close all;
