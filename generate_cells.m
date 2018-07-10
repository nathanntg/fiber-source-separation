function cells = generate_cells(varargin)
%GENERATE_CELLS Generate neural population
%   Given a volume, cell density and efficacy, generates uniformly
%   distributed neruons throughout the region,

% Size of the brain tissue volume being simulated
% from Bottjer el al 1985
% Area X is 2.262 mm^3
% assuming sphere, that works out to a radius of 814.3
% use slightly smaller X and Y dimensions
volume = [1200; 1200; 400]; % microns

% Density of cells to be simulated
% from Walton et al 2012, based on excitatory cells in HVC
% 275,000 per mm^3
cell_density = 0.000275; % cells per micron cubed

% Virus efficacy (adjusts cell density)
viral_efficacy = 1;

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

end
