function x = generate_inputs(number, p, dur, waveform, amplitude)
%GENERATE_INPUTS Summary of this function goes here
%   Detailed explanation goes here

% generate spike trains
x = binornd(1, p, number, dur);
% x = double(rand(number, duration) > p); <-- faster

% scale by amplitude
if isa(amplitude, 'function_handle')
    a = amplitude(numel(x));
    x = x .* reshape(a, size(x));
elseif number == length(amplitude)
    a = amplitude(:) * ones(1, dur);
    x = x .* a;
elseif 2 == length(amplitude)
    a = unifrnd(amplitude(1), amplitude(2), size(x, 1), size(x, 2));
    x = x .* a;
elseif 1 == length(amplitude)
    x = x * amplitude;
else
    error('Invalid value for amplitude.');
end

% convolve with waveform
x = conv2(x, waveform);

% trim (convolve adds extra values)
x = x(:, 1:dur);

end

