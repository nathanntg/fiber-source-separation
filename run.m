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
