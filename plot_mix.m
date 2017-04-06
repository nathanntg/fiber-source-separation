function b = plot_mix(mx, mode, display, split)
%PLOT_MIX Summary of this function goes here
%   Detailed explanation goes here

% defaults
if ~exist('mode', 'var')
    mode = 'most-total';
end
if ~exist('split', 'var')
    split = 7;
end
if ~exist('display', 'var')
    display = 50;
end

% counts
%fibers = size(mx, 1);
%cells = size(mx, 2);

% sort
switch mode
    case 'most-total'
        [~, idx] = sort(sum(mx, 1), 'descend');
    case 'most-single'
        [~, idx] = sort(max(mx, [], 1), 'descend');
    otherwise
        error('Invalid mode: "%s".' , mode);
end

% build b
b = zeros(display, split);
for j = 1:display
    i = idx(j);
    sc = sort(mx(:, i), 'descend');
    b(j, 1:(split - 1)) = sc(1:(split - 1));
    b(j, split) = sum(sc(split:end));
end

% display
if nargout == 0
    bar(b, 'stacked');
    ylabel('Normalized Fluence');
    xlim([0 display + 1]);
end

end

