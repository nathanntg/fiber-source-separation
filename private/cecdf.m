function [F, x] = cecdf(samples, bins)
%CEDF Continous Empirical Cumulative Distribution Function
%   This takes samples from a continuous distribution and forms an
%   empirical CDF with the specified number of bins (defaulting to 100). It
%   returns a CDF (F) and the associated x values (x). Note that if the
%   data has been discretized or rounded at all, bins should not exceed the
%   number of unique values in the sample.

% check that it exists
if ~exist('bins', 'var')
    bins = 100;
end

% get lower and upper bound
val_min = min(samples);
val_max = max(samples);
inc = (val_max - val_min) / bins;

% cnt, for normalization
cnt = max(size(samples));

% x values
x = val_min:inc:val_max;
F = zeros(size(x, 2), size(x, 1));

i = 1;
for j = x
    F(i) = sum(samples <= j) / cnt;
    i = i + 1;
end

% pad with 0, like ecdf does
x = [val_min; transpose(x)];
F = [0; F];

end

