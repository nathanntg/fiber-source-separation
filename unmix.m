function [M, X] = unmix(Y, waveform)
%UNMIX Summary of this function goes here
%   Detailed explanation goes here

% dimensions
waveform_len = length(waveform);
components = size(Y, 1);
timesteps = size(Y, 2);
timesteps_plus_padding = timesteps + waveform_len - 1;

% smoothing W matrix
W = zeros(timesteps_plus_padding, timesteps);
for i1 = 1:timesteps
    i2 = i1 + waveform_len - 1;
    W(i1:i2, i1) = waveform(end:-1:1);
end

% unsmooth
Y_deconvolved = Y / W(waveform_len:end, :);

% run NNMF for initial values
iterations = 1000;
opt = statset('MaxIter', iterations);
[M0, X0] = nnmf(Y_deconvolved, components, 'options', opt);

% pad X
X0 = [zeros(components, waveform_len - 1) X0];

% save as temporary file
save('temp.mat', '-v7', 'Y', 'M0', 'X0', 'waveform');

% call unmixing command
[status, cmd_out] = system('/Users/nathan/ve/unmix/bin/python unmix/unmix.py');
if status ~= 0
    fprintf('%s\n', cmd_out);
    error('Command failed.');
end

% load
res = load('temp.mat', 'M', 'Xw', 'costs', 'cost');
delete('temp.mat');

if nargout == 0
    % plot figure
    figure;
    semilogy(res.costs);
    title('Unmix Change');
    xlabel('Iteration');
    ylabel('M change; L2 norm');
end

% store results
M = res.M;
X = res.Xw;

end
