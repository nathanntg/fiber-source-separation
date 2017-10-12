%% SCENARIO 1: input / output
explore({}, 'number_of_inputs', 10:5:50, 'number_of_outputs', 10:5:75);
print(gcf, 'input-output.png', '-dpng', '-r300');
close all;

%% SCENARIO 2: input / noise
explore({'number_of_outputs', 25}, 'number_of_inputs', 10:5:50, 'input_noise', 0.1:0.15:0.90);
print(gcf, 'input-noise.png', '-dpng', '-r300');
close all;

%% SCENARIO 3: input / duration
explore({'number_of_outputs', 25}, 'number_of_inputs', 10:5:50, 'duration', 100:100:2000);
print(gcf, 'input-duration.png', '-dpng', '-r300');
close all;

%% SCENARIO 4: noise / duration
explore({'number_of_outputs', 25}, 'input_noise', 0.1:0.15:0.90, 'duration', 500:1500:10000);
print(gcf, 'noise-duration.png', '-dpng', '-r300');
close all;

%% SCENARIO 5: noise / amplitude
explore({'number_of_inputs', 25, 'number_of_outputs', 25}, 'input_noise', 0:0.25:2, 'amplitude', [1 2; 2 3; 4 5; 5 6; 7 8]);
print(gcf, 'noise-amplitude.png', '-dpng', '-r300');
close all;

%% SCENARIO 6: smoothing
explore({'number_of_outputs', 25}, 'number_of_inputs', 10:1:25, 'smooth_mixing', generate_smoothing(1, 1:9));
print(gcf, 'smoothing-avg-row.png', '-dpng', '-r300');
close all;

explore({'number_of_outputs', 25}, 'number_of_inputs', 10:1:25, 'smooth_mixing', generate_smoothing(1:9, 1));
print(gcf, 'smoothing-avg-column.png', '-dpng', '-r300');
close all;

explore({'number_of_outputs', 25}, 'number_of_inputs', 10:1:25, 'smooth_mixing', generate_smoothing(1:9, 1:9));
print(gcf, 'smoothing-avg-square.png', '-dpng', '-r300');
close all;

explore({'number_of_outputs', 25}, 'number_of_inputs', 10:1:25, 'smooth_mixing', generate_smoothing(1:9, 1:9, @(n, m) fspecial('gaussian', [n m])));
print(gcf, 'smoothing-gaus-square.png', '-dpng', '-r300');
close all;

explore({'number_of_outputs', 25}, 'number_of_inputs', 10:1:25, 'smooth_mixing', [{[]} generate_smoothing(3:2:9, 3:2:9, @(n, m) fspecial('disk', (n - 1) / 2))]);
print(gcf, 'smoothing-disk-square.png', '-dpng', '-r300');
close all;

%% SCENARIO 7: realistic output/duration
explore({'realistic', true}, 'number_of_outputs', 100:100:1200, 'duration', 1000:1000:15000, 3);
print(gcf, 'input-duration.png', '-dpng', '-r300');
close all;

%% SECNARIO 8: realistic ICA mode
profile = sp_model('sensitivity-profile/fiber-exc.mat');
profile = sp_3d_to_2d(profile); % symmetric, way faster
explore({'realistic', true, 'number_of_outputs', 250, 'duration', 2500, 'profile', profile}, 'g', {'pow3', 'tanh', 'gauss', 'skew'});
print(gcf, 'ica-g.png', '-dpng', '-r300');
close all;
