%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 4);
set(0, 'DefaultAxesFontSize', 24);

%% load
% excitation
profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
profile_exc = sp_3d_to_2d(profile_exc); % symmetric, way faster

% fluorescence 
profile_fluor = sp_model('sensitivity-profile/fiber-fluor.mat');
profile_fluor = sp_3d_to_2d(profile_fluor); % symmetric, way faster

% other settings
duration = 5000;
lm = 100;

%% open figure
figure;

all_scores = {};

%% no noise, 50 fibers
fibers = 50; noise = 0;
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 1);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('50 fibers; no noise');

%% no noise, 100 fibers
fibers = 100; noise = 0;
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 3);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('100 fibers; no noise');

%% no noise, 500 fibers
fibers = 500; noise = 0;
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 5);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('500 fibers; no noise');

%% no noise, 1000 fibers
fibers = 1000; noise = 0;
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 7);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('1000 fibers; no noise');

%% noise, 50 fibers
fibers = 50; noise = @(n) normrnd(0, 0.1, 1, n);
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 2);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('50 fibers; noise (\sigma = 0.1)');

%% noise, 100 fibers
fibers = 100; noise = @(n) normrnd(0, 0.1, 1, n);
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 4);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('100 fibers; noise (\sigma = 0.1)');

%% noise, 500 fibers
fibers = 500; noise = @(n) normrnd(0, 0.1, 1, n);
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 6);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('500 fibers; noise (\sigma = 0.1)');

%% noise, 1000 fibers
fibers = 1000; noise = @(n) normrnd(0, 0.1, 1, n);
scores = simulate_source_separation('mode', 'profile-rt', 'figures', false, ...
    'number_of_outputs', fibers, 'output_noise', noise, ...
    'duration', duration, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);
all_scores{end + 1} = scores;
subplot(4, 2, 8);
plot(sort(scores, 'descend'));
xlim([0 lm]); ylim([0 1]);
ylabel('r^2');
title('1000 fibers; noise (\sigma = 0.1)');
