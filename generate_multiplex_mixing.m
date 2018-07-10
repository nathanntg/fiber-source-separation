function m = generate_multiplex_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, mp_number, mp_all, varargin)
%GENERATE_MULTIPLEX_MIXING Multiplex mixing matrices
%   Given structured illumination of a subset of the fibers, you will end
%   up with distinct mixing matrices that may aid in signal reconstruction.

%% DEFAULT PARAMETERS
opt_fiber = [];
opt_neuron = [];
as_single = false;

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

% defaults
if ~exist('mp_all', 'var') || isempty(mp_all)
    mp_all = false;
end

%% GENERATE MIXING MATRIX

% number of fibers
number_of_outputs = size(fibers, 2);

if as_single
    m = single([]);
else
    m = [];
end

% if include all or need to prune cells?
if (mp_all && mp_number > 1) || ~isempty(opt_neuron)
    % generate regular mixing matrix
    new_m = generate_realistic_rt_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, 'figures', false, 'stats', false);
    
    % prune mixing columns and list of cells
    if ~isempty(opt_neuron)
        idx = any(new_m > opt_neuron, 1);
        cells = cells(:, idx);
        new_m = new_m(:, idx);
    end
    
    % include all?
    if mp_all && mp_number > 1
        % prune fibers
        if ~isempty(opt_fiber)
            new_m = new_m(any(new_m > opt_fiber, 2), :);
        end
        
        % append
        if as_single
            m = [m; single(new_m)];
        else
            m = [m; new_m];
        end
    end
end

% for each multiplexing scheme...
for i = 1:mp_number
    % generate limited fluoresence mixing matrix
    new_m = generate_realistic_rt_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, 'fibers_exc', i:mp_number:number_of_outputs, 'figures', false, 'stats', false);
    
    % prune fibers
    if ~isempty(opt_fiber)
        new_m = new_m(any(new_m > opt_fiber, 2), :);
    end
    
    % appned
    if as_single
        m = [m; single(new_m)];
    else
        m = [m; new_m];
    end
end

% re-prune neurons?
if ~isempty(opt_neuron)
    m = m(:, any(m > opt_neuron, 1));
end

end

