function ret = generate_inputs(number, freq, mps, sps, duration, waveform, offset, amplitude)
%GENERATE_INPUTS Summary of this function goes here
%   Detailed explanation goes here

% generate inputs in blocks to reduce memory usage
block_size = 10000;

% pad start
pad_start = numel(waveform);

% prepare output
ret = zeros(number, round(duration * sps));

% probability
p = freq / mps; % convert to frequency per model sample

% model duration
smpl = round(mps / sps); % subsample (approximation)
sample = 1 + ((0:(round(duration * sps) - 1)) * smpl);
dur_model = sample(end);

for i_start = 1:block_size:number
    i_end = min(i_start + block_size - 1, number);
    i_cnt = 1 + i_end - i_start;

    % generate spike trains
    x = binornd(1, p, i_cnt, pad_start + dur_model);
    % x = double(rand(number, duration) > p); <-- faster

    % scale by amplitude
    if isa(amplitude, 'function_handle')
        a = amplitude(numel(x));
        x = x .* reshape(a, size(x));
    elseif number == numel(amplitude)
        a = amplitude(i_start:i_end) * ones(1, pad_start + dur);
        x = x .* a;
    elseif 2 == numel(amplitude)
        a = unifrnd(amplitude(1), amplitude(2), size(x, 1), size(x, 2));
        x = x .* a;
    elseif isscalar(amplitude)
        x = x * amplitude;
    else
        error('Invalid value for amplitude.');
    end

    % convolve with waveform (use waveform as a digital filter)
    x = filter(waveform, 1, x, [], 2);

    % trim (remove start padding)
    x = x(:, (1 + pad_start):end);

    % add constant offset to each neuron
    if isa(offset, 'function_handle')
        o = offset(i_cnt);
        o = reshape(o, [], 1); % ensure column vector
    elseif number == numel(offset)
        o = reshape(offset(i_start:i_end), [], 1); % ensure column vector
    elseif 2 == numel(offset)
        o = normrnd(offset(1), offset(2), i_cnt, 1);
    elseif isscalar(offset)
        o = offset .* ones(i_cnt, 1);
    elseif isempty(offset)
        o = 0;
    else
        error('Invalid value for offset.');
    end

    % actually add offset
    x = bsxfun(@plus, x, o);
    
    ret(i_start:i_end, :) = x(:, sample);
end

end

