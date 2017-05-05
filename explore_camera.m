%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 3);
set(0, 'DefaultAxesFontSize', 16);

% store random number seed
old_rng = rng;

%% load
% excitation
profile_exc = sp_model('sensitivity-profile/camera-exc.mat');
profile_exc = sp_3d_to_2d(profile_exc); % symmetric, way faster

%% open plots
% plot excitation
im = [profile_exc.volume(end:-1:2, :); profile_exc.volume];
y = [-1 * profile_exc.r(end:-1:2) profile_exc.r];
figure;
imagesc(profile_exc.z, y, log10(im), [-4 0]);
title('Camera excitation profile');
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');

%% plot distributions for single fibers
rng(0);
[m_exc, ~] = generate_camera_mixing(1, 1, profile_exc, 'figures', false);
s = sort(m_exc, 'descend');
figure;
bar(s(1:150)); xlim([0 150]);
xlabel('Neuron'); ylabel('Normalized Fluence'); title('Excitation - 1 pixel');

%% average over multiple single fibers
iter = 50;
number = 150;
ss = zeros(iter, number);
for i = 1:iter
    [m_exc, ~] = generate_camera_mixing(1, 1, profile_exc, 'figures', false, 'stats', false);
    s = sort(m_exc, 'descend');
    ss(i, :) = s(1:number);
end
figure;
bar(mean(ss)); xlim([0 number + 1]);
xlabel('Neuron'); ylabel('Normalized Fluence'); title('Excitation - average per pixel');

%% plot distributions for multiple fibers
rng(0);
width = 20; height = 20;
[m_exc, ~] = generate_camera_mixing(width, height, profile_exc, 'figures', false);

% nice plot
figure;
plot_mix(m_exc);
xlabel('Neuron');
title(sprintf('Excitation - %d pixels', number));

% nice plot
figure;
plot_mix(m_exc, 'most-total', 500);
xlabel('Neuron');
title(sprintf('Excitation - %d pixels', number));

% nice plot 2
figure;
plot_mix(m_exc, 'most-single');
xlabel('Neuron');
title(sprintf('Excitation - %d pixels', number));

% nice plot
figure;
plot_mix(m_exc', 'most-single', 100);
xlabel('Pixel');
title('Excitation');

% plot kurtosis
figure;
plot(sort(kurtosis(m_exc, 0, 1), 'descend')); title('Kurtosis');
figure;
plot(sort(skewness(m_exc, 0, 1), 'descend')); title('Skewness');


%% explore number
areas = [50 125 200 275 350 425 500]; % SD
numbers = 100:100:1000;
well = 0.1;
seen = zeros(length(numbers), length(areas));
seen_well = zeros(length(numbers), length(areas));
seen_well_multiple = zeros(length(numbers), length(areas));
seen_well_aggregate = zeros(length(numbers), length(areas));
condition = zeros(length(numbers), length(areas));
for i = 1:length(numbers)
    for j = 1:length(areas)
        rng(0);
        [m, ~] = generate_realistic_mixing(numbers(i), profile_exc, 'fiber_distribution', [areas(j) 0 0; 0 areas(j) 0; 0 0 15], 'figures', false, 'stats', false);

        seen(i, j) = size(m, 2);
        seen_well(i, j) = sum(max(m, [], 1) > well);
        seen_well_multiple(i, j) = sum(sum(m > well, 1) > 1);
        seen_well_aggregate(i, j) = sum(sum(m, 1) > well);
        condition(i, j) = cond(m);
    end
end

%figure; plot(numbers, seen); title('Number of seen');
figure;
plot(numbers, seen_well);
title(sprintf('Cells seen well by one fiber (normalized fluence > %d%%)', well * 100));
h = legend(num2str(areas'), 'Location', 'bestoutside');
title(h, 'SD of Splay ({\sigma})');

figure;
plot(numbers, seen_well_multiple);
title(sprintf('Cells seen well by 2+ fibers (normalized fluence > %d%%)', well * 100));
h = legend(num2str(areas'), 'Location', 'bestoutside');
title(h, 'SD of Splay ({\sigma})');

figure;
plot(numbers, seen_well_aggregate);
title(sprintf('Cells seen well in aggregate (total normalized fluence > %d%%)', well * 100));
xlabel('Number of fibers');
h = legend(num2str(areas'), 'Location', 'bestoutside');
title(h, 'SD of Splay ({\sigma})');

figure;
plot(numbers, condition);
title('Condition number');
h = legend(num2str(areas'), 'Location', 'bestoutside');
title(h, 'SD of Splay ({\sigma})');

%% restore
rng(old_rng);
