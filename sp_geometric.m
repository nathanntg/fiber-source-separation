function profile = sp_geometric(varargin)
%SP_GEOMETRIC Summary of this function goes here
%   Detailed explanation goes here

% fiber width
fiber_diameter = 8; % microns
fiber_half_angle_of_acceptance = 16.2; % degrees

% core diameter
core_diameter = 2 * sqrt(0.4 * (fiber_diameter / 2) ^ 2);

% profile depth
profile_resolution = 0.5;
profile_depth = 500; % microns

%% LOAD PARAMETERS
nparams=length(varargin);
if 0 < mod(nparams, 2)
	error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

%% MAKE PROFILE

% slope of light cone
cone_radius = core_diameter / 2;
cone_slope = tan(deg2rad(fiber_half_angle_of_acceptance));

% depths
max_radius = cone_radius + profile_depth * cone_slope;
x = 0:profile_resolution:max_radius; x = [-1 * x(end:-1:2) x];
y = 0:profile_resolution:max_radius; y = [-1 * y(end:-1:2) y];
z = 0:profile_resolution:profile_depth;

% start with empty volume
volume = zeros(length(x), length(y), length(z));

% make mesh grid
[xx, yy] = meshgrid(x, y);
rr = xx .^ 2 + yy .^ 2;

% calculate radiuses at each depth
radius2 = (cone_radius + z .* cone_slope) .^ 2;

% fill volume
for i = 1:length(z)
    idx = rr <= radius2(i);
    strength = (cone_radius ^ 2) ./ radius2(i);
    v = volume(:, :, i);
    v(idx) = strength;
    volume(:, :, i) = v;
end

% build profile structure
profile = struct('x', x, 'y', y, 'z', z, 'volume', volume, 'd', 3);

end
