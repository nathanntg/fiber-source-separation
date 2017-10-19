function [m, exc_cells, exc_fibers] = generate_realistic_rt_mixing(fibers, fiber_angles, cells, profile_exc, profile_emi, varargin)

% works in 3D space where X and Y are perpindicular to the fibers,
% and Z is parallel to the fibers

% which fibers are used for excitaiton and which fibers are used for
% emissions
fibers_exc = []; % defaults to all fibers
fibers_emi = []; % defaults to all fibers

% fiber width
fiber_diameter = 8;

% show figures
figures = true;
stats = true;

%% LOAD PARAMETERS
nparams = length(varargin);
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

num_fibers = size(fibers, 2);
num_cells = size(cells, 2);
if isempty(fibers_exc)
    fibers_exc = 1:num_fibers;
end
if isempty(fibers_emi)
    fibers_emi = 1:num_fibers;
end

% make scatter histogram
if figures
    figure;
    scatterhist(fibers(1, :), fibers(2, :));
    title('Distribution of fibers');
end

%% APPLY EXCITATION PROFILE
exc_fibers = zeros(1, num_fibers); % how much light is at tip of each fiber (used for normalization)
exc_cells = zeros(1, num_cells); % how much light reaches each cell
for i = fibers_exc
    % get relative fiber positions
    fibers_rel = bsxfun(@minus, fibers, fibers(:, i));
    
    % get relative cell positions
    cells_rel = bsxfun(@minus, cells, fibers(:, i));
    
    % apply rotation
    theta = fiber_angles(:, i); % angle
    r = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))]; % x
    r = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))] * r; % y
    fibers_rel = r * fibers_rel; % rotate fibers
    cells_rel = r * cells_rel; % rotate cells
    
    % apply profile to fibers
    cur = apply_profile(profile_exc, fibers_rel);
    exc_fibers = exc_fibers + cur;
    
    % apply profile to cells
    cur = apply_profile(profile_exc, cells_rel);
    exc_cells = exc_cells + cur;
end

%% APPLY EMISSION PROFILE
m = zeros(num_fibers, num_cells);
for i = fibers_emi
    % get relative cell positions
    cells_rel = bsxfun(@minus, cells, fibers(:, i));
    
    % apply rotation
    theta = fiber_angles(:, i); % angle
    r = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))]; % x
    r = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))] * r; % y
    cells_rel = r * cells_rel; % rotate
    
    % apply profile
    cur = apply_profile(profile_emi, cells_rel);
    
    % add to mixing matrix
    m(i, :) = cur;
end

% total = excitation * emission
normalize = max(exc_fibers); % brighest point (at tip of fibers)
m = bsxfun(@times, m, exc_cells ./ normalize);

%% SUMMARY
if stats
    well = 0.1 * normalize;
    fprintf('Number of cells seen: %d\n', size(m, 2));
    fprintf('Number of cells seen well: %d\n', sum(max(m, [], 1) > well));
    fprintf('Cells seen by multiple fibers: %d\n', sum(sum(m > 0, 1) > 1));
    fprintf('Cells seen well by multiple fibers: %d\n', sum(sum(m > well, 1) > 1));
    %fprintf('Condition number: %f\n', cond(m));
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
