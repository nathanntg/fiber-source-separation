function [S, M, W] = unmix_nnica(Y, waveform)
%UNMIX_NNICA Use de-convolution and NN-ICA to perform source separation
%   If the fluorescent waveform is provided, an deconvolution filter is
%   first applied (as this can be applied to the already mixed signals) and
%   then NN-ICA is used to whiten the data. Once complete, the signal is
%   reconvolved with the waveform.

unfilter = ~isempty(waveform);

% unfilter
if unfilter
    waveform_len = length(waveform);
    delta = zeros(1, waveform_len * 3);
    delta(1) = 1;
    r = deconv(delta, waveform);
    Y = filter(r, 1, Y, [], 2);
end

% unmix
[S, M, W] = nnica(Y);

% update M (non-negative mixture)
rows_with_neg_mix = find(any(M < 0, 2));
if ~isempty(rows_with_neg_mix)
    % non-negative
    S = max(S, 0);
    
    % only bother with rows that have a negative mixture
    for i = 1:length(rows_with_neg_mix)
        row = rows_with_neg_mix(i);
        M(row, :) = lsqnonneg(S', Y(row, :)');
    end
    W = inv(M);
end

% smooth again
if unfilter
    S = filter(waveform, 1, S, [], 2);
end

end
