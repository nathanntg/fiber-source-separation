function [m, num_inputs] = generate_realistic_mixing(num_fibers, varargin)

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
fiber_width = 8;
fiber_half_angle_of_acceptance = 17.5;

% show figures
figures = true;

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


%% GENERATE SPACE

% number of cells
num_cells = round(prod(volume) * cell_density * viral_efficacy);

% cells are uniformly located throughout the space
cells = bsxfun(@times, rand(3, num_cells), volume);

% fibers are dsitributed with a multivariate normal distribution
% assume targeted to centers
fibers = mvnrnd(volume ./ 2, fiber_distribution, num_fibers)';
fiber_angles = deg2rad(mvnrnd([0, 0], angle_distribution, num_fibers))';

% slope of light cone
cone_radius = fiber_width / 2;
cone_slope = tan(deg2rad(fiber_half_angle_of_acceptance));

% make scatter histogram
if figures
    figure;
    scatterhist(fibers(1, :), fibers(2, :));
    title('Distribution of fibers');
end

%% FIGURE OUT CONE OF VIEW
m = zeros(num_fibers, num_cells);
for i = 1:num_fibers
    % get relative cell positions
    cells_rel = bsxfun(@minus, cells, fibers(:, i));
    
    % apply rotation
    theta = fiber_angles(:, i); % angle
    r = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))]; % x
    r = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))] * r; % y
    cells_rel = r * cells_rel; % rotate
    
    % calculate light cone radius
    radius2 = (cone_radius + cells_rel(3, :) * cone_slope) .^ 2;
    
    % in radius
    idx = (cells_rel(1, :) .^ 2 + cells_rel(2, :) .^ 2) <= radius2;
    idx = idx & cells_rel(3, :) >= 0;
    
    % constant strength
    strength = (cone_radius ^ 2) ./ radius2;
    
    % add to mixing matrix
    m(i, :) = idx .* strength;
end

% remove unused cells
m = m(:, any(m > 0, 1));
num_inputs = size(m, 2);

%% SUMMARY
    
well = 0.1;
fprintf('Number of cells seen: %d\n', num_inputs);
fprintf('Number of cells seen well: %d\n', sum(max(m, [], 1) > well));
fprintf('Cells seen by multiple fibers: %d\n', sum(sum(m > 0, 1) > 1));
fprintf('Cells seen well by multiple fibers: %d\n', sum(sum(m > well, 1) > 1));

if figures
    figure;
    imagesc(m);
    title('Mixing matrix');
    
    figure;
    hist(m(:), 100);
    title('Distribution of weights');
end

end