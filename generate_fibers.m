function [fibers, fiber_angles] = generate_fibers(num_fibers, varargin)
%GENERATE_FIBERS Generate fiber distribution
%   Given a volume, center position and splay parameters (normal
%   distribution standard deviation), generates a splayed bundle of fibers
%   consistent with previously measured histology.

% Size of the brain tissue volume being simulated
% from Bottjer el al 1985
% Area X is 2.262 mm^3
% assuming sphere, that works out to a radius of 814.3
% use slightly smaller X and Y dimensions
volume = [1200; 1200; 400]; % microns

% Center of fiber bundle implant.
position = [0.5; 0.5; 0.25];

% Distribution of fibers is assumed to be normal
% LNY 63 has bivariate 49.7-95.4 micron standard deviation (assuming
% independent), use slightly larger for more fibers
fiber_distribution = [125 0 0; 0 125 0; 0 0 15]; % microns

% angle of fibers, in degrees
angle_distribution = 0; % [5 0; 0 5];

% fiber width
fiber_diameter = 8;

% fiber distribution
distribution = 'mv';

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

if isscalar(fiber_distribution)
    fiber_distribution = [1 0 0; 0 1 0; 0 0 0.12] * fiber_distribution;
end
if isscalar(angle_distribution)
    angle_distribution = [1 0; 0 1] * angle_distribution;
end

%% GENERATE

% sanity check fiber count
fibers_1sd = num_fibers * 0.682; % number of fibers in 1 standard deviation
fiber_area_1sd = fibers_1sd * pi * (fiber_diameter / 2) ^ 2; % area of all the fibers
area_1sd = pi * fiber_distribution(1, 1) * fiber_distribution(2, 2); % only handles independent bivariate
circle_packing = 0.9069;
if fiber_area_1sd > (circle_packing * area_1sd)
    warning('Based on circle packing and the fiber density, the maximum number of fibers is %d.', round((circle_packing * area_1sd) / (pi * (fiber_diameter / 2) ^ 2) / 0.682));
end

% fiber position
if isempty(position)
    position = volume ./ 2;
else
    if isscalar(position)
        position = ones(size(volume));
    end
    
    % allow fraction
    a = position > 0 & position < 1;
    position(a) = position(a) .* volume(a);
end

if strcmp(distribution, 'uniform')
    % uniform distribution is useful for some situations where we want to
    % directly measure the effect of varying distribution parameters
    fibers = [...
        unifrnd(position(1) - fiber_distribution(1, 1) / 2, position(1) + fiber_distribution(1, 1) / 2, 1, num_fibers), ...
        unifrnd(position(2) - fiber_distribution(2, 2) / 2, position(2) + fiber_distribution(2, 2) / 2, 1, num_fibers), ...
        unifrnd(position(3) - fiber_distribution(3, 3) / 2, position(3) + fiber_distribution(3, 3) / 2, 1, num_fibers) ...
        ];
    fiber_angles = deg2rad([...
        unifrnd(0 - angle_distribution(1, 1) / 2, 0 + angle_distribution(1, 1) / 2, 1, num_fibers), ...
        unifrnd(0 - angle_distribution(2, 2) / 2, 0 + angle_distribution(2, 2) / 2, 1, num_fibers) ...
        ]);
else
    % fibers are dsitributed with a multivariate normal distribution
    % assume targeted to centers (mvnrnd takes VARIANCE)
    fibers = mvnrnd(position, fiber_distribution .^ 2, num_fibers)';
    fiber_angles = deg2rad(mvnrnd([0, 0], angle_distribution .^ 2, num_fibers))';
end

% check bounds
if any(min(fibers, [], 2) < 0) || any(max(fibers, [], 2) > volume)
    warning('Fibers exceed modeled volume.');
end

end

