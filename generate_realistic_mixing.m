function [m, num_inputs] = generate_realistic_mixing(num_fibers, profile, varargin)

% works in 3D space where X and Y are perpindicular to the fibers,
% and Z is parallel to the fibers

% Size of the brain tissue volume being simulated
% from Bottjer el al 1985
% Area X is 2.262 mm^3
% assuming sphere, that works out to a radius of 814.3
% use slightly smaller X and Y dimensions
volume = [1200; 1200; 400];

% Density of cells to be simulated
% from Walton et al 2012, based on excitatory cells in HVC
% 275,000 per mm^3
cell_density = 0.000275; % cells per micron cubed

% Virus efficacy (adjusts cell density)
viral_efficacy = 1;

% Distribution of fibers is assumed to be normal
% LNY 63 has bivariate 49.7-95.4 micron standard deviation (assuming
% independent), use slightly larger for more fibers
fiber_distribution = [125 0 0; 0 125 0; 0 0 15];

% angle of fibers, in degrees
angle_distribution = [5 0; 0 5];

% fiber width
fiber_diameter = 8;

% show figures
figures = true;
stats = true;

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

%% GENERATE SPACE

% number of cells
num_cells = round(prod(volume) * cell_density * viral_efficacy);

% cells are uniformly located throughout the space
cells = bsxfun(@times, rand(3, num_cells), volume);

% sanity check fiber count
fibers_1sd = num_fibers * 0.682; % number of fibers in 1 standard deviation
fiber_area_1sd = fibers_1sd * pi * (fiber_diameter / 2) ^ 2; % area of all the fibers
area_1sd = pi * fiber_distribution(1, 1) * fiber_distribution(2, 2); % only handles independent bivariate
circle_packing = 0.9069;
if fiber_area_1sd > (circle_packing * area_1sd)
    warning('Based on circle packing and the fiber density, the maximum number of fibers is %d.', round((circle_packing * area_1sd) / (pi * (fiber_diameter / 2) ^ 2)));
end

% fibers are dsitributed with a multivariate normal distribution
% assume targeted to centers
fibers = mvnrnd(volume ./ 2, fiber_distribution, num_fibers)';
fiber_angles = deg2rad(mvnrnd([0, 0], angle_distribution, num_fibers))';


% make scatter histogram
if figures
    figure;
    scatterhist(fibers(1, :), fibers(2, :));
    title('Distribution of fibers');
end

%% APPLY PROFILE
m = zeros(num_fibers, num_cells);
for i = 1:num_fibers
    % get relative cell positions
    cells_rel = bsxfun(@minus, cells, fibers(:, i));
    
    % apply rotation
    theta = fiber_angles(:, i); % angle
    r = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))]; % x
    r = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))] * r; % y
    cells_rel = r * cells_rel; % rotate
    
    % calculate strength from sensitivity profile
    cur = interp3(profile.x, profile.y, profile.z, profile.volume, cells_rel(1, :), cells_rel(2, :), cells_rel(3, :));
    
    % outside
    idx = cells_rel(1, :) < profile.x(1) | cells_rel(1, :) > profile.x(end) | ...
        cells_rel(2, :) < profile.y(1) | cells_rel(2, :) > profile.y(end) | ...
        cells_rel(3, :) < profile.z(1) | cells_rel(3, :) > profile.z(end);
    
    % clear outside
    cur(idx) = 0;
    
    % add to mixing matrix
    m(i, :) = cur;
end

% remove unused cells
m = m(:, any(m > 0, 1));
num_inputs = size(m, 2);

%% SUMMARY

if stats
    well = 0.1;
    fprintf('Number of cells seen: %d\n', num_inputs);
    fprintf('Number of cells seen well: %d\n', sum(max(m, [], 1) > well));
    fprintf('Cells seen by multiple fibers: %d\n', sum(sum(m > 0, 1) > 1));
    fprintf('Cells seen well by multiple fibers: %d\n', sum(sum(m > 0, 1) > 1 & max(m, [], 1) > well));
    fprintf('Condition number: %f\n', cond(m));
end

if figures
    figure;
    imagesc(m);
    title('Mixing matrix');
    
    figure;
    hist(m(:), 100);
    title('Distribution of weights');
    
    % core diameter
    core_diameter = 2 * sqrt(0.4 * (fiber_diameter / 2) ^ 2);
    
    % slope of light cone
    fiber_half_angle_of_acceptance = 16.2; % degrees
    cone_radius = core_diameter / 2;
    cone_slope = tan(deg2rad(fiber_half_angle_of_acceptance));
    
    figure;
    c = jet(256);
    scatter3(cells_rel(1, :), cells_rel(2, :), cells_rel(3, :), 50, c(round(1 + mat2gray(cur) .* 210 + (cur > 0) .* 24), :), 'filled');
    [x, y, z] = cylinder([cone_radius; cone_radius + 60 * cone_slope]);
    hold on;
    o = surf(x, y, z * 60);
    alpha(o, 0.1);
    hold off;
    set(gca,'ZDir','reverse');
    xlim([-50 50]);
    ylim([-50 50]);
    zlim([-25 75]);
    xlabel('{\mu}m'); ylabel('{\mu}m'); zlabel('{\mu}m');
end

end
