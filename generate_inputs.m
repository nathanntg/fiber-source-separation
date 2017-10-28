function x = generate_inputs(number, p, dur, waveform, offset, amplitude)
%GENERATE_INPUTS Summary of this function goes here
%   Detailed explanation goes here

% generate spike trains
x = binornd(1, p, number, dur);
% x = double(rand(number, duration) > p); <-- faster

% scale by amplitude
if isa(amplitude, 'function_handle')
    a = amplitude(numel(x));
    x = x .* reshape(a, size(x));
elseif number == numel(amplitude)
    a = amplitude(:) * ones(1, dur);
    x = x .* a;
elseif 2 == numel(amplitude)
    a = unifrnd(amplitude(1), amplitude(2), size(x, 1), size(x, 2));
    x = x .* a;
elseif isscalar(amplitude)
    x = x * amplitude;
else
    error('Invalid value for amplitude.');
end

% convolve with waveform
x = conv2(x, waveform);

% trim (convolve adds extra values)
x = x(:, 1:dur);

% add constant offset to each neuron
if isa(offset, 'function_handle')
    o = offset(number);
    o = reshape(o, [], 1); % ensure column vector
elseif number == numel(offset)
    o = reshape(offset, [], 1); % ensure column vector
elseif 2 == numel(offset)
    o = normrnd(offset(1), offset(2), number, 1);
elseif isscalar(offset)
    o = offset .* ones(number, 1);
elseif isempty(offset)
    o = 0;
else
    error('Invalid value for offset.');
end

% actually add offset
x = bsxfun(@plus, x, o);

end

