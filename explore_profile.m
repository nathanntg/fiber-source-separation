%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 3);
set(0, 'DefaultAxesFontSize', 16);

% store random number seed
old_rng = rng;

%% load
% excitation
profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');
profile_exc = sp_3d_to_2d(profile_exc); % symmetric, way faster

% fluorescence 
profile_fluor = sp_model('sensitivity-profile/fiber-fluor.mat');
profile_fluor = sp_3d_to_2d(profile_fluor); % symmetric, way faster

%% open plots
% plot excitation
im = [profile_exc.volume(end:-1:2, :); profile_exc.volume];
y = [-1 * profile_exc.r(end:-1:2) profile_exc.r];
figure;
imagesc(profile_exc.z, y, log10(im), [-4 0]);
title('Excitation profile');
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');

% plot fluoresence
im = [profile_fluor.volume(end:-1:2, :); profile_fluor.volume];
y = [-1 * profile_fluor.r(end:-1:2) profile_fluor.r];
figure;
imagesc(profile_fluor.z, y, log10(im), [-4 0]);
title('Fluoresence profile');
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');

%% plot distributions for single fibers
rng(0);
[fibers, fiber_angles] = generate_fibers(1);
cells = generate_cells();
[m_exc, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false);
s = sort(m_exc, 'descend');
figure;
bar(s(1:150)); xlim([0 150]);
xlabel('Neuron'); ylabel('Normalized Fluence'); title('Excitation - 1 fiber');

%% average over multiple single fibers
iter = 50;
number = 100;
ss = zeros(iter, number);
a = [];
for i = 1:iter
    [fibers, fiber_angles] = generate_fibers(1);
    cells = generate_cells();
    [m_exc, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false, 'stats', false);
    s = sort(m_exc, 'descend');
    ss(i, :) = s(1:number);
    a = [a sum(s > 0.1)];
end
figure;
bar(mean(ss)); xlim([0 number + 1]);
xlabel('Neuron'); ylabel('Normalized Fluence'); title('Excitation - average per fiber');

%% plot distributions for multiple fibers
rng(0);
number = 1000;
[fibers, fiber_angles] = generate_fibers(number);
cells = generate_cells();
[m_exc, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false);

% nice plot
figure;
plot_mix(m_exc);
xlabel('Neuron');
title(sprintf('Excitation - %d fibers', number));

% nice plot
figure;
plot_mix(m_exc, 'most-total', 500);
xlabel('Neuron');
title(sprintf('Excitation - %d fibers', number));

% nice plot 2
figure;
plot_mix(m_exc, 'most-single');
xlabel('Neuron');
title(sprintf('Excitation - %d fibers', number));

% nice plot
figure;
plot_mix(m_exc', 'most-single', 100);
xlabel('Fiber');
title('Excitation');

% plot kurtosis
figure;
plot(sort(kurtosis(m_exc, 0, 1), 'descend')); title('Kurtosis');
figure;
plot(sort(skewness(m_exc, 0, 1), 'descend')); title('Skewness');

%% signal to background
signal = max(m_exc, [], 2);
background = sum(m_exc, 2) - signal;

% depths - scatter plot
[~, clearest_cell] = max(m_exc, [], 2);
s_signal = max(m_exc, [], 1);
figure;
depths = cells(3, unique(clearest_cell));
scatter(depths, s_signal(unique(clearest_cell)));
xlabel('Depth ({\mu})'); ylabel('Visibility');
title('Most Visible Neurons');

% depths - scatter plot
figure;
scatter(cells(3, :), s_signal);
xlabel('Depth ({\mu})'); ylabel('Visibility');
title('All Neurons');

% depth vs ratio - scatter plot
depths = [];
s_cbr = signal ./ background;
s2_cbr = [];
for j = unique(clearest_cell(:)')
    depths = [depths cells(3, j)];
    s2_cbr = [s2_cbr max(s_cbr(clearest_cell == j))];
end
figure;
scatter(depths, s2_cbr);
xlabel('Depth ({\mu})'); ylabel('Contrast to Background Ratio');

% distribution
figure;
[f, x] = cecdf(s_cbr);
plot(x, f);
xlabel('Contrast to Background Ratio (signal / background)');
%legend('Fiber', 'Microscope', 'Location', 'SouthEast');

%% explore unique cells
number = [100 200 500 1000 2000 5000];
brightest = [];
for j = number
    rng(0);
    [fibers, fiber_angles] = generate_fibers(j);
    cells = generate_cells();
    [m_exc, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false, 'stats', false);
    [~, clearest_cell] = max(m_exc, [], 2);
    brightest = [brightest length(unique(clearest_cell))];
end
figure;
scatter(number, brightest);
xlabel('Number of Fibers'); ylabel('Number of Bright Unique Neurons');

% adjust distribution
number = [100 200 500 1000 2000 5000];
brightest = [];
for j = number
    rng(0);
    [fibers, fiber_angles] = generate_fibers(j, 'fiber_distribution', [sqrt(j) * 5 0 0; 0 sqrt(j) * 5 0; 0 0 15]);
    cells = generate_cells();
    [m_exc, ~, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false, 'stats', false);
    [~, clearest_cell] = max(m_exc, [], 2);
    brightest = [brightest length(unique(clearest_cell))];
end
figure;
scatter(number, brightest);
xlabel('Number of Fibers'); ylabel('Number of Bright Unique Neurons');

%figure;
%plot(number, brightest, number1, brightest1, number2, brightest2);
%xlabel('Number of Fibers / Pixels'); ylabel('Number of Bright Unique Neurons');
%legend('Fibers', 'Fibers (constant splay)', 'Microscope', 'Location', 'SouthEast');

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
        [fibers, fiber_angles] = generate_fibers(numbers(i), 'fiber_distribution', [areas(j) 0 0; 0 areas(j) 0; 0 0 15]);
        cells = generate_cells();
        [m, ~] = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, 'figures', false, 'stats', false);

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
title(h, 'SD of Splay ({\mu}m)');
xlabel('Number of fibers');
ylabel('Number of neurons');

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
