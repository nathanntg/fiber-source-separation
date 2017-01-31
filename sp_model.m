function profile = sp_model(fn)
%SP_MODEL Summary of this function goes here
%   Detailed explanation goes here

%% READ FILE
a = load(fn);
I = a.I;
config = a.config;

%% ENCODE
% units for axes
x = linspace(...
    double(config.voxel_min(1) - 1) * config.voxel_size(1) * 1000,...
    double(config.voxel_max(1) - 1) * config.voxel_size(1) * 1000,...
    double(config.voxel_count(1)));
y = linspace(...
    double(config.voxel_min(2) - 1) * config.voxel_size(2) * 1000,...
    double(config.voxel_max(2) - 1) * config.voxel_size(2) * 1000,...
    double(config.voxel_count(2)));
z = linspace(...
    double(config.voxel_min(3) - 1) * config.voxel_size(3) * 1000,...
    double(config.voxel_max(3) - 1) * config.voxel_size(3) * 1000,...
    double(config.voxel_count(3)));

% adjust positions relative to source
x = x - config.initial_position(1) * 1000;
y = y - config.initial_position(2) * 1000;
z = z - config.initial_position(3) * 1000;

% make volume
volume = I;
volume(volume < 0) = 0;
volume = volume ./ max(volume(:));

% build profile structure
profile = struct('x', x, 'y', y, 'z', z, 'volume', volume);

end
