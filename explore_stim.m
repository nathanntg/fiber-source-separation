%% setup
% sizing
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 4);
set(0, 'DefaultAxesFontSize', 24);

%% prep variables
% excitation
profile_exc = sp_model('sensitivity-profile/fiber-exc.mat');

% calculate voxel size
voxel_size = median(diff(profile_exc.x)) * median(diff(profile_exc.y)) * median(diff(profile_exc.z));

% stimulation threshold
light_out = 0.0544  * 52; % uW out of the fiber.
light_out_per_mm = light_out * 15.83; % mW/mm^2; based on: "(1 µW) / (.4 * (? * 8 µm / 2) ^ 2) in mW/mm^2"
threshold = 5; % mW / mm^2

% cell density
density = 0.000275; % cells per micron cubed

%% single fiber
% above stimulation threshold
stim = profile_exc.volume * light_out_per_mm  > threshold;

% convert to cubic microns
stim_volume = sum(stim(:)) * voxel_size;
fprintf('%d um^3\n', stim_volume);
fprintf('Cells: %.1f\n', stim_volume * density);

% make slice
idx_zero = find(profile_exc.y == 0); % zero Y value
slice = squeeze(stim(idx_zero, :, :));

imagesc(profile_exc.z, profile_exc.x, slice)
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
title(sprintf('Stimulation threshold (%d mW / mm^2)', threshold));

%% setup multiple fibers fibers
number_of_fibers = 500;
volume = [1200; 1200; 1200];

% make fibers and cells
[fibers, fiber_angles] = generate_fibers(number_of_fibers, ...
    'position', [0.5; 0.5; 100], 'volume', volume);
cells = generate_cells('volume', volume, 'cell_density', density);
m = generate_realistic_mixing(fibers, fiber_angles, cells, profile_exc, ...
    'stats', false, 'figures', false);

% scale m
m_scaled = m .* light_out_per_mm;

% potentially drop cells that would never be stimulated?
drop = sum(m .* light_out_per_mm, 1) < threshold;
m_scaled(:, drop) = [];

%% calculate multiple fibers
% calculate number activated by single fibers
above_thresh = sum(m_scaled > threshold, 2);

% debug
fprintf('Average cells per fiber: %.1f\n', mean(above_thresh));

% full population
ca_individual = sum(any(m_scaled > threshold, 1));
ca_cumulative = sum(sum(m_scaled, 1) > threshold);
fprintf('Activated individually: %d\n', ca_individual);
fprintf('Activated cumulatively: %d\n', ca_cumulative);

% try various numbers of fibers
subsets = [5 10 25 50 100];
iterations = 20000;
output = zeros(iterations, length(subsets), 2);
for i = 1:length(subsets)
    num = subsets(i);
    for j = 1:iterations
        idx = randperm(number_of_fibers, num);
        
        % calculate cells activated ("ca")
        
        % total given individual
        ca_individual = sum(any(m_scaled(idx, :) > threshold, 1));
        
        % total given cumulative power
        ca_cumulative = sum(sum(m_scaled(idx, :), 1) > threshold);
        
        output(j, i, :) = [ca_individual, ca_cumulative];
    end
end

% activated
fprintf('Number Activated:\n');
disp(mean(output(:, :, 2) - output(:, :, 1)));

% percentage activated
fprintf('Number Activated:\n');
disp(mean((output(:, :, 2) - output(:, :, 1)) ./ output(:, :, 1)));

