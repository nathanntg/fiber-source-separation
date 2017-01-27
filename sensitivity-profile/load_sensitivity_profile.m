function [I, IoOut] = load_sensitivity_profile(fn, config)

% autoload configuration
if ~exist('config', 'var') || isempty(config)
    [d, f, ~] = fileparts(fn);
    config = load_configuration(fullfile(d, [f '.inp']));
end

% is double precision
s = dir(fn);
bytes = s(1).bytes;
threshold = prod([config.voxel_count config.angular_steps([2 1])]) * 6;
if bytes > threshold
    units = 'float64';
else
    units = 'float32';
end

% load second 2pt
fid = fopen(fn);
Io = fread(fid, units);
fclose(fid);

% separate 5d matrix and photon exit data
IoOut = Io((end-1):end);
Io = reshape(Io(1:end-2), [config.voxel_count config.angular_steps([2 1])]);

% normalize the 2pt
lstIn = find(Io > 0);
lstOut = find(Io < 0);
Io(lstOut) = Io(lstOut) / config.photons;
IoOut(2) = IoOut(2) / config.photons;
nPhotonDet = -(sum(Io(lstOut)) + IoOut(2));
k = (1 - nPhotonDet) / (config.tissue_absorption(1) * prod(config.voxel_size) * (sum(Io(lstIn)) + IoOut(1)));
Io(lstIn) = Io(lstIn) * k;
IoOut(1) = IoOut(1) * k;

% integrate over angles
I = sum(sum(Io, 5), 4);

end
