function output = explore(params, param_a_name, param_a_values, param_b_name, param_b_values, iterations, detailed, metric)

% params
if isempty(params)
    params = {};
end
params = [params {'figures', false}];

% defaults
if ~exist('iterations', 'var') || isempty(iterations)
    iterations = 5;
end
if nargin < 4 || isempty(param_b_name)
    is2D = false;
    param_b_values = 0;
else
    is2D = true;
end
if ~exist('detailed', 'var') || isempty(detailed)
    detailed = false;
end

if detailed
    output = cell(length(param_a_values), length(param_b_values));
else
    output = zeros(length(param_a_values), length(param_b_values), iterations);
end

if ~exist('metric', 'var') || isempty(metric)
    metric = 0.7;
else
    if detailed
        error('Metric only used if not detailed mode.');
    end
    if ischar(metric) && ~strcmp(metric, 'mean') && ~strcmp(metric, 'isi') && ~strcmp(metric, 'auc')
        error('The only supported metrics are a threshold (scalar value) or "mean", "isi" or "auc".');
    end
end

% progress
cnt = 0;
max_cnt = numel(output);

% explore space
for idx_a = 1:length(param_a_values)
    val_a = param_a_values(idx_a);
    for idx_b = 1:length(param_b_values)
        val_b = param_b_values(idx_b);
        
        % add parameters and values
        if is2D
            p = [params param_a_name val_a param_b_name val_b];
        else
            p = [params param_a_name val_a];
        end
        
        % make holder for output values
        values = zeros(1, iterations);
        
        % run number of iterations
        for i = 1:iterations
            [scores, isi, auc] = simulate_source_separation(p{:});
            if detailed
                output{idx_a, idx_b, i} = scores;
            elseif ischar(metric) && strcmp(metric, 'mean') % mean
                values(i) = mean(scores);
            elseif ischar(metric) && strcmp(metric, 'isi') % isi
                values(i) = isi;
            elseif ischar(metric) && strcmp(metric, 'auc') % auc
                values(i) = auc;
            else
                values(i) = sum(scores > metric);
            end
        end
        
        % store output
        if ~detailed
            output(idx_a, idx_b, :) = values;
        end
        
        % print progress
        cnt = cnt + iterations;
        fprintf('Progress: %d%%\n', round(100 * cnt / max_cnt));
    end
end

% plot figure
if nargout < 1 && ~detailed
    figure;
    if is2D
        x = param_b_values;
        if iscell(x)
            x = 1:length(param_b_values);
        elseif 1 < min(size(x))
            x = x(:, 1);
        end
        y = param_a_values;
        if iscell(y)
            y = 1:length(param_a_values);
        elseif 1 < min(size(y))
            y = y(:, 1);
        end
        if ~ischar(metric) && isscalar(metric)
            number_of_outputs = 50;
            if strcmp(param_a_name, 'number_of_outputs')
                number_of_outputs = max(param_a_values);
            elseif strcmp(param_b_name, 'number_of_outputs')
                number_of_outputs = max(param_b_values);
            else
                for i = 1:2:length(params)
                    if strcmp(params{i}, 'number_of_outputs')
                        number_of_outputs = params{i + 1};
                    end
                end
            end
            im = mean(output, 3);
            imagesc(x, y, im, [0 min(number_of_outputs, max(im(:)))]);
            title(sprintf('# of separated signals (r^2 \\geq %.1f)', metric));
        else
            imagesc(x, y, mean(output, 3));
            title('Accuracy of source separation');
        end
        axis xy;
        xlabel(strrep(param_b_name, '_' ,' '));
        ylabel(strrep(param_a_name, '_' ,' '));
        colorbar;
    else
        bar(mean(output, 3));
        set(gca, 'XTickLabel', param_a_values);
        xlabel(strrep(param_a_name, '_' ,' '));
        if ~ischar(metric) && isscalar(metric)
            ylabel(sprintf('# of separated signals (r^2 \\geq %.1f)', metric));
        else
            ylabel('Accuracy of source separation');
        end
        
        figure;
        boxplot(squeeze(output)');
        set(gca, 'XTickLabel', param_a_values);
        xlabel(strrep(param_a_name, '_' ,' '));
        if ~ischar(metric) && isscalar(metric)
            ylabel(sprintf('# of separated signals (r^2 \\geq %.1f)', metric));
        else
            ylabel('Accuracy of source separation');
        end
    end
end

end
