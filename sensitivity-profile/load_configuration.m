function config = load_configuration(fn)
%LOAD_CONFIGURATION Summary of this function goes here
%   Detailed explanation goes here

config = struct();

% load input file
fid = fopen(fn);

% load number of photons
c = textscan(fid, '%d\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected number of photons');
end
config.photons = double(c{1});

% load random seed
c = textscan(fid, '%d\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected random seed');
end
config.seed = double(c{1});

% load initial position
c = textscan(fid, '%f %f %f\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected initial position');
end
config.initial_position = c{1};

% load inital direction
c = textscan(fid, '%f %f %f\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected initial direction');
end
config.initial_direction = c{1};

% load time steps
c = textscan(fid, '%f %f %f\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected time steps');
end
config.time_range = cat(2, c{1, 1:2});
config.time_step = c{1, 3};

% load angles
c = textscan(fid, '%d %d\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected angles');
end
config.angular_steps = double(cat(2, c{1, 1:2}));

% load tissue file
c = textscan(fid, '%s\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected tissue file');
end
config.tissue_file = c{1};

% load dimensions
c = textscan(fid, '%f %d %d %d\n', 3);
if size(c, 1) ~= 1 || length(c{1}) ~= 3
    fclose(fid);
    error('Unexpected volume dimensions');
end
config.voxel_size = cat(2, c{1})';
config.volume_step = c{2}';
config.voxel_min = c{3}';
config.voxel_max = c{4}';
config.voxel_count = 1 + config.voxel_max - config.voxel_min;

% number of tissue types
c = textscan(fid, '%d\n', 1);
if size(c, 1) ~= 1
    fclose(fid);
    error('Unexpected number of tissue types');
end
tissue_types = c{1};

% load tissue types
c = textscan(fid, '%f %f %f %f\n', tissue_types);
if size(c, 1) ~= tissue_types
    fclose(fid);
    error('Unexpected tissue type data');
end
config.tissue_scattering = c{1}';
config.tissue_anisotropy = c{2}';
config.tissue_absorption = c{3}';
config.tissue_refractive = c{4}';

% close file
fclose(fid);

end
