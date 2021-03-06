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

%% ONE OFF
simulate_source_separation('mode', 'multiplex-profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor);

%% AUC
simulate_source_separation('mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 100, 'number_of_outputs', 50, 'auc_threshold', 0.2:0.1:0.8);
print(gcf, 'auc.png', '-dpng', '-r150');
close all;

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
% effect of simulation duration
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'duration', 50:50:300);
print(gcf, 'duration.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: NOISE
% effect of noise
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250}, 'output_noise', 0:0.05:0.2);
print(gcf, 'noise.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: SPIKE FREQUENCY
% effect of spike frequency
explore({'mode', 'profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250}, 'spike_frequency', 0.02:0.4:2.2);
print(gcf, 'spike_frequency.png', '-dpng', '-r150');
close all;

%% SCRATCH PAD: MULTIPLEX
explore({'mode', 'multiplex-profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250, 'number_of_outputs', 1000}, 'mp_number', [1 5 10 50 100 500]);
print(gcf, 'multiplex.png', '-dpng', '-r150');

explore({'mode', 'multiplex-profile-rt', 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'duration', 250, 'number_of_outputs', 1000}, 'mp_number', [1 5 10 50 100 500], 'g', 'unmix');
print(gcf, 'multiplex-unmix.png', '-dpng', '-r150');
%close all;

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

%% SECNARIO 10: samples per second
explore({'mode', 'profile-rt', 'number_of_outputs', 100, 'duration', 250, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor}, 'sps', 5:5:25);
print(gcf, 'sps.png', '-dpng', '-r300');
close all;

%% SECNARIO 11: multiplex
explore({'mode', 'multiplex-profile-rt', 'duration', 250, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'sps', 10, 'mp_all', false}, 'mp_number', 1:4, 'number_of_outputs', 50:50:250);
print(gcf, 'multiplex.png', '-dpng', '-r300');
close all;
