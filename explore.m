function output = explore(params, param_a_name, param_a_values, param_b_name, param_b_values, iterations)

% params
if isempty(params)
    params = {};
end
params = [params {'figures', false}];

% defaults
if nargin < 6
    iterations = 10;
end
if nargin < 4 || isempty(param_b_name)
    is2D = false;
    param_b_values = 0;
else
    is2D = true;
end

output = zeros(length(param_a_values), length(param_b_values));

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
            scores = simulate_source_separation(p{:});
            values(i) = mean(scores);
        end
        
        % store output
        output(idx_a, idx_b) = mean(values);
        
        % print progress
        cnt = cnt + 1;
        fprintf('Progress: %d%%\n', round(100 * cnt / max_cnt));
    end
end

% plot figure
if nargout < 1
    figure;
    if is2D
        x = param_b_values;
        if 1 < min(size(x))
            x = x(:, 1);
        end
        y = param_a_values;
        if 1 < min(size(y))
            y = y(:, 1);
        end
        imagesc(x, y, output);
        axis xy;
        title('Accuracy of source separation');
        xlabel(strrep(param_b_name, '_' ,' '));
        ylabel(strrep(param_a_name, '_' ,' '));
        colorbar;
    else
        bar(param_a_values, output);
        xlabel(strrep(param_a_name, '_' ,' '));
        ylabel('Average correlation');
    end
end

end
