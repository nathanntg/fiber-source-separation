function [scores, isi, auc] = simulate_source_separation(varargin)

%#ok<*UNRCH>

% mode
sss_mode = ''; % special modes for generating figures: not a great approach

% simulation settings
duration = 50; % seconds
mps = 100; % model time steps per second
sps = 20; % samples per second
spike_frequency = 0.4; % Hz (spikes per seconds)
number_of_inputs = 25; % n: sources or neurons
number_of_outputs = 50; % m: number of distringuishable fibers
% mixing matrix will m x n size

% waveform (inputs are spike trains convolved with a waveform)
waveform = 'gcamp6f'; % the waveform measured from each spike (vector or string)
offset = [0.1 0.01]; % waveform offset (scalar = constant, pair of parameters [mu, sigma], or callback)
amplitude = [1 2]; % either single value (all the same amplitude), or range (uniform), or callback

% noise
input_noise = 0; % either single value (0) or callback that generates noise
input_noise_type = 'scale';
output_noise = @(n) normrnd(0, 0.1, 1, n); % either single value (0) or callback that generates noise
output_noise_type = 'scale';

figures = true;

% correlations
smooth_mixing = []; % row: nearby inputs influence output; column: nearby outputs influence each other

% mixing
mode = 'basic'; % basic, profile, profile-rt
profile = []; % used for profile
profile_exc = []; % used for profile-rt
profile_fluor = []; % used for profile-rt

% params fibers
params_fibers = {};
params_cells = {};

% multiplex
mp_number = 2;
mp_all = false;

% ICA
g = 'skew';

% threshold (optimization, ignore weak cells)
opt_threshold = 1e-7; % eps

% auc threshold
auc_threshold = 0.7;

%% LOAD PARAMETERS
nparams = length(varargin);
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
switch mode
    case 'multiplex-profile-rt'
        % default profile (geometric)
        if isempty(profile_exc)
            profile_exc = sp_geometric();
        end
        if isempty(profile_fluor)
            profile_fluor = profile_exc; % use same one
        end
        
        % prepare profile and other parameters
        [fibers, fiber_angles] = generate_fibers(number_of_outputs, params_fibers{:});
        cells = generate_cells(params_cells{:});
        
        % generate mixing matrix
        if mp_number > 100
            warning('Using heavily optimized mixing matrix.');
            m = generate_multiplex_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, mp_number, mp_all, 'opt_fiber', 1e-3, 'opt_neuron', 1e-4);
        else
            m = generate_multiplex_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, mp_number, mp_all, 'opt_fiber', 1e-3, 'opt_neuron', 1e-4);
        end
        
        % remove unused cells
        m = m(:, any(m > opt_threshold, 1));
        number_of_inputs = size(m, 2);
        
    case 'profile-rt'
        % default profile (geometric)
        if isempty(profile_exc)
            profile_exc = sp_geometric();
        end
        if isempty(profile_fluor)
            profile_fluor = profile_exc; % use same one
        end
        
        % prepare profile and other parameters
        [fibers, fiber_angles] = generate_fibers(number_of_outputs, params_fibers{:});
        cells = generate_cells(params_cells{:});
        m = generate_realistic_rt_mixing(fibers, fiber_angles, cells, profile_exc, profile_fluor, 'figures', false, 'stats', false);

        % remove unused cells
        m = m(:, any(m > opt_threshold, 1));
        number_of_inputs = size(m, 2);
        
    case 'profile'
        % default profile (geometric)
        if isempty(profile)
            profile = sp_geometric();
        end
        
        % prepare profile and other parameters
        [fibers, fiber_angles] = generate_fibers(number_of_outputs, params_fibers{:});
        cells = generate_cells(params_cells{:});
        m = generate_realistic_mixing(fibers, fiber_angles, cells, profile, 'figures', false, 'stats', false);

        % remove unused cells
        m = m(:, any(m > opt_threshold, 1));
        number_of_inputs = size(m, 2);
        
    case 'basic'
        m = generate_mixing_matrix(number_of_inputs, number_of_outputs);

        % remove unused cells
        m = m(:, any(m > opt_threshold, 1));
        number_of_inputs = size(m, 2);
        
    otherwise
        error('Invalid mode: %s.', mode);
end

if ~isempty(smooth_mixing)
    m = filter2(smooth_mixing, m);
end

%% GENERATE INPUT
waveform_v = get_waveform(waveform, mps); % wave form in model samples
s = generate_inputs(number_of_inputs, spike_frequency, mps, sps, duration, waveform_v, offset, amplitude);
duration_smp = size(s, 2); % duration in readout samples
s_noisy = add_noise(s, input_noise, input_noise_type);

% show example
if figures
    figure;
    subplot(3, 1, 1);
    plot((1:200) ./ sps, s_noisy(1, 1:200));
    title('Example source traces');
    xlabel('Time');
    ylabel('Trace (s_1)');
    subplot(3, 1, 2);
    plot((1:200) ./ sps, s_noisy(2, 1:200));
    xlabel('Time');
    ylabel('Trace (s_2)');
    subplot(3, 1, 3);
    plot((1:200) ./ sps, s_noisy(3, 1:200));
    xlabel('Time');
    ylabel('Trace (s_3)');
end

%% GENERATE OUTPUTS
x = m * s_noisy;
x_noisy = add_noise(x, output_noise, output_noise_type);
if figures
    figure;
    subplot(3, 1, 1);
    plot((1:200) ./ sps, x(1, 1:200)', (1:200) ./ sps, x_noisy(1, 1:200)');
    title('Example output traces');
    xlabel('Time');
    ylabel('Trace (x_1)');
    subplot(3, 1, 2);
    plot((1:200) ./ sps, x(2, 1:200)', (1:200) ./ sps, x_noisy(2, 1:200)');
    xlabel('Time');
    ylabel('Trace (x_2)');
    subplot(3, 1, 3);
    plot((1:200) ./ sps, x(3, 1:200)', (1:200) ./ sps, x_noisy(3, 1:200)');
    xlabel('Time');
    ylabel('Trace (x_3)');
end

%% PERFORM SOURCE SEPARATION
if strcmp(g, 'unmix')
    [m_hat, s_hat] = unmix(x_noisy, get_waveform(waveform, sps));
    w_hat = inv(m_hat);
elseif strcmp(g, 'nnmf')
    [m_hat, s_hat] = unmix_nnmf(x_noisy, get_waveform(waveform, sps));
    w_hat = inv(m_hat);
elseif strcmp(g, 'none')
    s_hat = x_noisy;
    m_hat = eye(number_of_outputs);
    w_hat = eye(number_of_outputs);
else
    % perform ICA
    if figures
        % s_hat, m_hat, w_hat
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'g', g, 'numOfIC', number_of_inputs);
    else
        % s_hat, m_hat, w_hat
        [s_hat, m_hat, w_hat] = fastica(x_noisy, 'verbose', 'off', 'displayMode', 'off', 'g', g, 'numOfIC', number_of_inputs);
    end

    % no convergence?
    if isempty(s_hat)
        % perform ICA
        if figures
            % s_hat, m_hat, w_hat
            [s_hat, m_hat, w_hat] = fastica(x_noisy, 'numOfIC', number_of_inputs);
        else
            % s_hat, m_hat, w_hat
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
end

%% ANALYZE RESULTS
% score outcome
rho = corr(s_noisy', s_hat');

% figure out best scores
[scores, idx] = max(abs(rho), [], 2);

% special mode for paper_model.m
if strcmp(sss_mode, 'traces')
    % sort
    [~, sorted_in] = sort(scores);
    sorted_out = idx(sorted_in);
    
    % brightest cells in mixing matrix
    [~, brightest_per_fiber] = sort(m_hat, 2, 'descend');
    best_scores = max(rho, [], 1);
    
    % find where the brightest two cells both have scores about 0.7
    potential_fibers = find(all(best_scores(brightest_per_fiber(:, 1:2)) > 0.7, 2));
    
    % brightest
    brightest_intensity_per_fiber = sort(m_hat(potential_fibers, :), 2, 'descend');
    [~, ratio_between_brightest] = min(brightest_intensity_per_fiber(:, 1) ./ brightest_intensity_per_fiber(:, 2));
    
    % fiber to plot
    fiber_to_plot = potential_fibers(ratio_between_brightest);
    
    % strongest neurons
    [~, separated_in_fiber] = sort(m_hat(fiber_to_plot, :), 'descend');
    [~, neurons_in_fiber] = max(rho(:, separated_in_fiber));
    
    clrs = lines(3);

    h = figure('Renderer', 'painters');
    nm = min(duration_smp, 300);
    t = (1:nm) ./ sps;
    
    subplot(3, 1, 1);
    plot(t, mat2gray(x_noisy(fiber_to_plot, 1:nm)), 'Color', clrs(2, :)); % fiber signal
    ylim([0 1]); yticks([]);
    xlabel('Time [s]'); ylabel('Intensity');
    title(sprintf('Fiber (y_{%d})', fiber_to_plot));
    legend('Fiber', 'Location', 'NorthWest');
    
    subplot(3, 1, 2);
    hold on;
    plot(t, 0.5 + mat2gray(s(neurons_in_fiber(1), 1:nm)), 'Color', clrs(1, :)); % original signal
    plot(t, 0 + mat2gray(s_hat(separated_in_fiber(1), 1:nm)), 'Color', clrs(3, :)); % separated signal
    hold off;
    ylim([0 1.5]); yticks([]);
    xlabel('Time [s]'); ylabel('Trace');
    title(sprintf('Neuron (x_{%d})', neurons_in_fiber(1)));
    legend('Neuron', 'Separated', 'Location', 'NorthWest');
    
    subplot(3, 1, 3);
    hold on;
    plot(t, 0.5 + mat2gray(s(neurons_in_fiber(2), 1:nm)), 'Color', clrs(1, :)); % original signal
    plot(t, 0 + mat2gray(s_hat(separated_in_fiber(2), 1:nm)), 'Color', clrs(3, :)); % separated signal
    hold off;
    ylim([0 1.5]); yticks([]);
    xlabel('Time [s]'); ylabel('Trace');
    title(sprintf('Neuron (x_{%d})', neurons_in_fiber(2)));
    
    h.Position(4) = h.Position(4) * 3;
    
    return;
end

if figures
    % sort
    [~, sorted_in] = sort(scores);
    sorted_out = idx(sorted_in);
    
    % plot best two
    figure;
    nm = min(duration_smp, 300);
    t = (1:nm) ./ sps;
    subplot(3, 3, 1);
    plot(t, s(sorted_in(end), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 2);
    plot(t, s_noisy(sorted_in(end), 1:nm));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 3);
    plot(t, s_hat(sorted_out(end), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 4);
    plot(t, s(sorted_in(end - 1), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 5);
    plot(t, s_noisy(sorted_in(end - 1), 1:nm));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 6);
    plot(t, s_hat(sorted_out(end - 1), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end - 1)));
    xlabel('Time'); ylabel('Trace');
    % plot worst
    subplot(3, 3, 7);
    plot(t, s(sorted_in(1), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 8);
    plot(t, s_noisy(sorted_in(1), 1:nm));
    title(sprintf('Source (s_{%d}) + Noise', sorted_in(1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 3, 9);
    plot(t, s_hat(sorted_out(1), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(1)));
    xlabel('Time'); ylabel('Trace');
    
    % plot best two
    figure;
    nm = min(duration_smp, 300);
    t = (1:nm) ./ sps;
    subplot(3, 2, 1);
    plot(t, s(sorted_in(end), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 2, 2);
    plot(t, s_hat(sorted_out(end), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 2, 3);
    plot(t, s(sorted_in(end - 1), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(end - 1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 2, 4);
    plot(t, s_hat(sorted_out(end - 1), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(end - 1)));
    xlabel('Time'); ylabel('Trace');
    % plot worst
    subplot(3, 2, 5);
    plot(t, s(sorted_in(1), 1:nm));
    title(sprintf('Source (s_{%d})', sorted_in(1)));
    xlabel('Time'); ylabel('Trace');
    subplot(3, 2, 6);
    plot(t, s_hat(sorted_out(1), 1:nm));
    title(sprintf('Separated Source (s_{%d})', sorted_out(1)));
    xlabel('Time'); ylabel('Trace');
    
    figure;
    plot(sort(scores, 'descend'));
    xlim([1 100]); ylim([0 1]);
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
        t = (1:duration_smp) ./ sps;
        
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
    
%     figure;
%     fpr = [];
%     tpr = [];
%     for r = 0:(length(sorted_out) - 1)
%         num = length(unique(sorted_out((end - r):end)));
%         fpr = [fpr (number_of_outputs - num) / number_of_outputs];
%         tpr = [tpr num / number_of_outputs];
%         if num == size(s_hat, 1)
%             break;
%         end
%     end
end

if nargout > 1
    % get best matches
    idx_in = [];
    for i = 1:size(s_hat, 1)
        [~, j] = max(rho(:, i));
        idx_in = [idx_in j];
    end
    
    % isi
    % global matrix (should be close to identity with a bunch of zeroish
    % rows, i guess)
    g = w_hat * m; % (:, idx_in);
    g_abs = abs(g);
    K = size(g, 1);
    isi = 0;
    for i = 1:size(g, 1)
        isi = isi + sum(g_abs(i, :)) / max(g_abs(i, :)) - 1;
    end
    for j = 1:size(g, 2)
        isi = isi + sum(g_abs(:, j)) / max(g_abs(:, j)) - 1;
    end
    isi = isi / (2 * K * (K - 1));
    
    % auc control
    if strcmp(sss_mode, 'auc_control')
        s = generate_inputs(number_of_inputs, spike_frequency, mps, sps, duration, waveform_v, offset, amplitude);
        sss_mode = 'auc';
    end
    
    % auc
    if strcmp(sss_mode, 'auc')
        % struct structure
        auc = struct('threshold', nan, 'number', nan, 'auc', nan, 'tpr', nan, 'fpr', nan);
        for i = 1:length(auc_threshold)
            rec = sum(scores > auc_threshold(i));
            [a, tpr, fpr] = calculate_auc(s, s_hat, scores, idx, auc_threshold(i));
            auc(i) = struct('threshold', auc_threshold(i), 'number', rec, 'auc', a, 'tpr', tpr, 'fpr', fpr);
        end
    elseif isscalar(auc_threshold)
        auc = calculate_auc(s, s_hat, scores, idx, auc_threshold);
    else
        auc = zeros(size(auc_threshold));
        for i = 1:length(auc_threshold)
            auc(i) = calculate_auc(s, s_hat, scores, idx, auc_threshold(i));
        end
    end
end

end
