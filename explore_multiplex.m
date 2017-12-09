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

%% explore density

% no multiplexing
explore({'mode', 'profile-rt', 'duration', 500, 'sps', 10, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'params_fibers', {'fiber_distribution', [50 0 0; 0 50 0; 0 0 5]}}, 'number_of_outputs', 50:50:200);
title('No multiplexing');

% yes multiplexing
explore({'mode', 'multiplex-profile-rt', 'duration', 500, 'sps', 10, 'profile_exc', profile_exc, 'profile_fluor', profile_fluor, 'mp_all', false, 'mp_number', 2, 'params_fibers', {'fiber_distribution', [50 0 0; 0 50 0; 0 0 5]}}, 'number_of_outputs', 50:50:200);
title('Multiplexing');
