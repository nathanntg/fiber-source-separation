function [M, X] = unmix_nnmf(Y, waveform)
%UNMIX_NNMF Summary of this function goes here
%   Detailed explanation goes here

unfilter = ~isempty(waveform);

% unfilter
if unfilter
    waveform_len = length(waveform);
    delta = zeros(1, waveform_len * 3);
    delta(1) = 1;
    r = deconv(delta, waveform);
    Y = filter(r, 1, Y, [], 2);
    Y = max(Y, 0);
end

% dimensions
components = size(Y, 1);

% run NNMF
iterations = 1000;
opt = statset('MaxIter', iterations);
[M, X] = nnmf(Y, components, 'options', opt);


% smooth again
if unfilter
    X = filter(waveform, 1, X, [], 2);
end

end

