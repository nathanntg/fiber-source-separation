% generate mixing matrix
[fibers, fiber_angles] = generate_fibers(100);
cells = generate_cells();
m = generate_realistic_mixing(fibers, fiber_angles, cells);

% total optical power to each cell
power_total = sum(m, 1);

% maxiumum optical power to each cell
power_max = max(m, [], 1);

% difference
power_diff = power_total - power_max;

% sort
[~, idx] = sort(power_diff, 'descend');

% bar graph
figure;
subplot(2, 2, 1);
bar([power_max(idx)' power_diff(idx)'], 'stacked');
ylabel('Irradiance'); xlim([1 length(power_diff)]);
subplot(2, 2, 3);
bar([power_max(idx)' power_diff(idx)'], 'stacked');
ylabel('Irradiance [log]'); xlim([1 length(power_diff)]);
set(gca, 'Yscale', 'log');

% line plot
subplot(2, 2, 2);
plot(1:length(power_total), power_total(idx), 1:length(power_max), power_max(idx));
ylabel('Irradiance'); xlim([1 length(power_diff)]);
subplot(2, 2, 4);
plot(1:length(power_total), power_total(idx), 1:length(power_max), power_max(idx));
ylabel('Irradiance [log]'); xlim([1 length(power_diff)]);
set(gca, 'Yscale', 'log');