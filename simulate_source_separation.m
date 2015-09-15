function [scores] = simulate_source_separation(varargin)

% simulation settings
duration = 1000; % number of samples
spike_probability = 0.02; % ~1 spike per 50 units of time
number_of_inputs = 25; % n: sources or neurons
number_of_outputs = 50; % m: number of distringuishable fibers
% mixing matrix will m x n size

% waveform (inputs are spike trains convolved with a waveform)
waveform = exp(0:-0.75:-4); % the waveform measured from each spike
amplitude = [1 2]; % either single value (all the same amplitude), or range (uniform), or callback

% noise
input_noise = @(n) normrnd(0, 0.25, 1, n); % either single value (0) or callback that generates noise
output_noise = 0; % either single value (0) or callback that generates noise

figures = true;

% correlations
% TODO: write me

%% LOAD PARAMETERS
nparams=length(varargin);
if 0 < mod(nparams, 2)
	error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

%% GENERATE INPUT
s = generate_inputs(number_of_inputs, spike_probability, duration, waveform, amplitude);
% TODO: add correlations to inputs
s_noisy = add_noise(s, input_noise);

% show example
if figures
    figure;
    subplot(3, 1, 1);
    plot(s_noisy(1, 1:200));
    title('Example source traces');
    xlabel('Time');
    ylabel('Trace (s_1)');
    subplot(3, 1, 2);
    plot(s_noisy(2, 1:200));
    xlabel('Time');
    ylabel('Trace (s_2)');
    subplot(3, 1, 3);
    plot(s_noisy(3, 1:200));
    xlabel('Time');
    ylabel('Trace (s_3)');
end

%% GENERATE MIXING MATRIX
m = generate_mixing_matrix(number_of_inputs, number_of_outputs);
% TODO: add correlations to mixing matrix

% show mixing matrix
if figures
    figure;
    imagesc(m);
    title('Mixing Matrix');
    xlabel('Inputs');
    ylabel('Outputs');
    colorbar;
end

%% GENERATE OUTPUTS
x = m * s_noisy;
x_noisy = add_noise(x, output_noise);
if figures
    figure;
    subplot(3, 1, 1);
    plot(x_noisy(1, 1:200)');
    title('Example output traces');
    xlabel('Time');
    ylabel('Trace (x_1)');
    subplot(3, 1, 2);
    plot(x_noisy(2, 1:200)');
    xlabel('Time');
    ylabel('Trace (x_2)');
    subplot(3, 1, 3);
    plot(x_noisy(3, 1:200)');
    xlabel('Time');
    ylabel('Trace (x_3)');
end

%% PERFORM SOURCE SEPARATION
% perform ICA
if figures
    [s_hat, m_hat, w_hat] = fastica(x_noisy, 'g', 'skew', 'numOfIC', number_of_inputs);
else
    [s_hat, m_hat, w_hat] = fastica(x_noisy, 'verbose', 'off', 'displayMode', 'off', 'g', 'skew', 'numOfIC', number_of_inputs);
end

% no convergence?
if isempty(s_hat)
    fprintf('WARNING: No convergence using g=skew.\n');
    
    % perform ICA
    if figures
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'numOfIC', number_of_inputs);
    else
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'verbose', 'off', 'displayMode', 'off', 'numOfIC', number_of_inputs);
    end
    
    % give up
    if isempty(s_hat)
        fprintf('WARNING: No convergence using g=pow3.\n');
        scores = 0;
        return;
    end
end

%% ANALYZE RESULTS
% score outcome
rho = corr(s_noisy', s_hat');
if figures
    figure;
    imagesc(rho);
    title('Correlation between sources and extracted sources');
    xlabel('Extracted source');
    ylabel('Source');
    colorbar;
end

% figure out best scores
[scores, idx] = max(rho, [], 2);
[~, sorted_in] = sort(scores);
sorted_out = idx(sorted_in);
if figures
    % plot best two
    subplot(2, 3, 1);
    plot(s(sorted_in(end), 1:200));
    title(sprintf('Source (s_{%d})', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(2, 3, 2);
    plot(s_noisy(sorted_in(end), 1:200));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(2, 3, 3);
    plot(s_hat(sorted_out(end), 1:200));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(2, 3, 4);
    plot(s(sorted_in(end - 1), 1:200));
    title(sprintf('Source (s_{%d})', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(2, 3, 5);
    plot(s_noisy(sorted_in(end - 1), 1:200));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(2, 3, 6);
    plot(s_hat(sorted_out(end - 1), 1:200));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end - 1)));
    xlabel('Time'); ylabel('Trace');
end

end