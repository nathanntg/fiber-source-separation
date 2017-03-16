function profile = sp_model(fn, threshold)
%SP_MODEL Summary of this function goes here
%   Detailed explanation goes here

%% DEFAULTS
if ~exist('threshold', 'var')
    threshold = []; % 1e-3;
end

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

% apply threshold
if ~isempty(threshold)
    above_threshold = volume > threshold;
    
    inc_x = any(any(above_threshold, 1), 3);
    inc_y = any(any(above_threshold, 2), 3);
    inc_z = any(any(above_threshold, 1), 2);
    
    x = x(inc_x);
    y = y(inc_y);
    z = z(inc_z);
    volume = volume(inc_y, inc_x, inc_z);
end

% build profile structure
profile = struct('x', x, 'y', y, 'z', z, 'volume', volume, 'd', 3);

end
