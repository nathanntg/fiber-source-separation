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
input_noise = @(n) normrnd(0, 0.125, 1, n); % either single value (0) or callback that generates noise
output_noise = 0; % either single value (0) or callback that generates noise

figures = true;

% correlations
smooth_input = [];
smooth_mixing = []; % row: nearby inputs influence output; column: nearby outputs influence each other

% mixing
realistic = false;
profile = [];

% ICA
g = 'skew';

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

%% GENERATE MIXING MATRIX
if realistic
    % prepare profile and other parameters
    [fibers, fiber_angles] = generate_fibers(number_of_outputs);
    cells = generate_cells();
    if isempty(profile)
        profile = sp_geometric();
    end
    
    m = generate_realistic_mixing(fibers, fiber_angles, cells, profile, 'figures', false, 'stats', false);
    number_of_inputs = size(m, 2);
    
    
    % remove unused cells
    m = m(:, any(m > 0, 1));
    number_of_inputs = size(m, 2);
else
    m = generate_mixing_matrix(number_of_inputs, number_of_outputs);
    
    % remove unused cells
    m = m(:, any(m > 0, 1));
end

if ~isempty(smooth_mixing)
    m = filter2(smooth_mixing, m);
end

% show mixing matrix
if figures
    figure;
    imagesc(m);
    title('Mixing Matrix');
    xlabel('Inputs');
    ylabel('Outputs');
    colorbar;
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
    [s_hat, m_hat, w_hat] = fastica(x_noisy, 'g', g, 'numOfIC', number_of_inputs);
else
    [s_hat, m_hat, w_hat] = fastica(x_noisy, 'verbose', 'off', 'displayMode', 'off', 'g', g, 'numOfIC', number_of_inputs);
end

% no convergence?
if isempty(s_hat)
    % perform ICA
    if figures
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'numOfIC', number_of_inputs);
    else
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'verbose', 'off', 'displayMode', 'off', 'numOfIC', number_of_inputs);
    end
    
    % give up
    if isempty(s_hat)
        fprintf('WARNING: No convergence using g=pow3. Giving up.\n');
        scores = 0;
        return;
    else
        fprintf('WARNING: No convergence using g=skew. Used g=pow3 instead.\n');
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
    figure;
    subplot(3, 3, 1);
    plot(s(sorted_in(end), 1:200));
    title(sprintf('Source (s_{%d})', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 2);
    plot(s_noisy(sorted_in(end), 1:200));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 3);
    plot(s_hat(sorted_out(end), 1:200));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 4);
    plot(s(sorted_in(end - 1), 1:200));
    title(sprintf('Source (s_{%d})', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 5);
    plot(s_noisy(sorted_in(end - 1), 1:200));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 6);
    plot(s_hat(sorted_out(end - 1), 1:200));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end - 1)));
    xlabel('Time'); ylabel('Trace');
    % plot worst
    subplot(3, 3, 7);
    plot(s(sorted_in(1), 1:200));
    title(sprintf('Source (s_{%d})', sorted_in(1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 8);
    plot(s_noisy(sorted_in(1), 1:200));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 9);
    plot(s_hat(sorted_out(1), 1:200));
    title(sprintf('Separated Source (s_{%d})', sorted_out(1)));
    xlabel('Time'); ylabel('Trace');
    
    figure;
    plot(sort(scores, 'descend'));
    xlim([1 length(scores)]); ylim([0 1]);
    title('Scores');
    ylabel('r');
    xlabel('Number of signals');
    
    % plot time courses
    if number_of_outputs <= 100
        % TODO: figure out better order
        
        % get best matches
        idx_in = [];
        idx_out = 1:size(s_hat, 1);
        for i = idx_out
            [~, j] = max(rho(:, i));
            idx_in = [idx_in j];
        end
        
        % colors
        clrs = lines(max(number_of_inputs, number_of_outputs));
        t = 1:duration;
        
        figure;
        
        % get traces
        ss = s_noisy(idx_in, :);
        ts = bsxfun(@minus, ss, min(ss, [], 2));
        ts = bsxfun(@rdivide, ts, max(ts, [], 2));
        [is, ~] = find(ts > 0.5 & cumsum(ts > 0.5, 2) == 1); % fist entry in each row
        % because of the way linear indexing works, `is` is in order
        is = is(end:-1:1);
        % actually plot
        subplot(1, 2, 1);
        hold on;
        xlim([t(1) t(end)]); ylim([0 size(ss, 1)]); xlabel('Time (s)');
        set(gca,'ytick',[]); title('Original Signals');
        for j = 1:size(ss, 1)
            i = is(j);
            %i = j;
            trace = ts(i, :);
            hold on; plot(t, j - 1 + trace, 'Color', clrs(j, :));
        end
        hold off;
        
        % get traces
        ss = s_hat(idx_out, :);
        ts = bsxfun(@minus, ss, min(ss, [], 2));
        ts = bsxfun(@rdivide, ts, max(ts, [], 2));
        %[is, ~] = find(ts > 0.5 & cumsum(ts > 0.5, 2) == 1); % fist entry in each row
        % because of the way linear indexing works, `is` is in order
        %is = is(end:-1:1);
        %actually plot
        subplot(1, 2, 2);
        xlim([t(1) t(end)]); ylim([0 size(ss, 1)]); xlabel('Time (s)');
        set(gca,'ytick',[]); title('Separated Signals');
        hold on;
        for j = 1:size(ss, 1)
            i = is(j);
            trace = ts(i, :);
            hold on; plot(t, j - 1 + trace, 'Color', clrs(j, :));
        end
        hold off;
    end
end

end