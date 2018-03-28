function [M, Xw] = unmix_nnmf(Y, waveform)
%UNMIX_NNMF Summary of this function goes here
%   Detailed explanation goes here

% dimensions
components = size(Y, 1);
timesteps = size(Y, 2);
waveform_len = length(waveform);

% smoothing W matrix
W = zeros(timesteps, timesteps);
for i2 = 1:timesteps
    i1 = max(i2 - waveform_len + 1, 1);
    l = 1 + i2 - i1;
    W(i1:i2, i2) = waveform(l:-1:1);
end

% unsmooth
Y = Y / W;

% run NNMF
iterations = 1000;
opt = statset('MaxIter', iterations);
[M, X] = nnmf(Y, components, 'options', opt);

% re-smooth
Xw = X * W;

end

