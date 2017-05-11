% load configuration
config = load_configuration('fiber-fluor.inp');

% load sensitivity profile
[I, IoOut] = load_sensitivity_profile('fiber-fluor.2pt', config);

% save for faster loading
save('fiber-fluor.mat', '-v7.3', 'I', 'config');

% get slice through source position
iFy = round(config.initial_position(2) * size(I, 2));

% units for axes
x = linspace(...
    double(config.voxel_min(1) - 1) * config.voxel_size(1),...
    double(config.voxel_max(1) - 1) * config.voxel_size(1),...
    double(config.voxel_count(1)));
z = linspace(...
    double(config.voxel_min(3) - 1) * config.voxel_size(3),...
    double(config.voxel_max(3) - 1) * config.voxel_size(3),...
    double(config.voxel_count(3)));

% make image
im = log10(squeeze(max(I(:, iFy, :), 1e-10)));

% show image
figure;
% convert units to microns
imagesc(z .* 1000, x .* 1000, im, [0 4]);
colormap(jet);
% add axis and labels
axis xy; xlabel('z [{\mu}]'); ylabel('x [{\mu}]');
