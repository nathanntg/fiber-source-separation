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

%% SCRATCH PAD: OFFSET
% effect of ofset
%explore({'mode', 'profile-rt', 'number_of_outputs', 25, 'duration', 100, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'offset', {0, 0.1, [0.1 0.01], 1, [1 1.1]});

%% SCRATCH PAD: WAVEFORM
% effect of waveform
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 500}, 'waveform', {'gcamp6s', 'gcamp6m', 'gcamp6f'});
print(gcf, 'waveform.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: OUTPUT NOISE TYPE
% effect of waveform
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250}, 'output_noise_type', {'add', 'scale'});
print(gcf, 'output_noise_type.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: DURATION
% effect of waveform
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'duration', 50:50:300);
print(gcf, 'duration.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: NOISE
% effect of waveform
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250}, 'output_noise', 0:0.05:0.2);
print(gcf, 'noise.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: P
% effect of waveform
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250}, 'spike_probability', 0.005:0.005:0.03);
print(gcf, 'spike_probability.png', '-dpng', '-r150');
close all;

%% SCENARIO 1: fibers / duration
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'number_of_outputs', 50:50:500, 'duration', 50:50:250);
print(gcf, 'fibers-duration.png', '-dpng', '-r300');
close all;

%% SCENARIO 2: fibers / noise
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'number_of_outputs', 50:50:300, 'output_noise', 0:0.05:0.2);
print(gcf, 'fibers-noise.png', '-dpng', '-r300');
close all;

%% SCENARIO 4: noise / duration
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'output_noise', 0:0.05:0.2, 'duration', 50:50:250);
print(gcf, 'noise-duration.png', '-dpng', '-r300');
close all;

%% SCENARIO 5: noise / amplitude
explore({'number_of_inputs', 25, 'number_of_outputs', 25}, 'input_noise', 0:0.25:2, 'amplitude', [1 2; 2 3; 4 5; 5 6; 7 8]);
print(gcf, 'noise-amplitude.png', '-dpng', '-r300');
close all;

%% SECNARIO 9: realistic round-trip, ICA mode
explore({'mode', 'profile-rt', 'number_of_outputs', 100, 'duration', 250, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'g', {'pow3', 'tanh', 'gauss', 'skew'});
print(gcf, 'ica-g.png', '-dpng', '-r300');
close all;
