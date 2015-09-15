function m_noisy = add_noise(m, noise)
%ADD_NOISE Adds noise to a signal. Takes a range of noise options.
%   Can take a function handle to compute any noise (function must take a
%   single parameter for number of noise elements to generate).
%   Alternatively, can take one or two parameters for normally distributed
%   noise (if two are provided, mu and sigma, otherwise just sigma and mu
%   defaults to zero). Can also take 0 for noise and will add no noise.

% implement different types of noise
if isa(noise, 'function_handle')
    % any function can be used to generate noise, must accept number of
    % elements
    n = noise(numel(m));
    m_noisy = m + reshape(n, size(m));
elseif 2 == length(noise)
    % 2 parameters is assumed to be mu and sigma for normally distributed
    % noise
    n = normrnd(noise(1), noise(2), size(m));
    m_noisy = m + n;
elseif 1 == length(noise)
    % 1 parameter is either constant or assumed to be sigma for zero-mean
    % normal noise
    if noise == 0
        m_noisy = m;
    else
        n = normrnd(0, noise, size(m));
        m_noisy = m + n;
    end
else
    error('Invalid value for noise.');
end

end

