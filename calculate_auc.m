function [auc, tpr, fpr, th] = calculate_auc(s, s_hat, scores, idx, r_threshold)
%CALCULATE_AUC Summary of this function goes here
%   Detailed explanation goes here

row_orig = find(scores > r_threshold);
row_hat = idx(scores > r_threshold);

% no matches?
if isempty(row_hat)
    auc = nan;
    tpr = [];
    fpr = [];
    th = nan;
    return;
end

orig = zeros(0, 'logical');
hat = [];

for i = 1:length(row_orig)
    % append original
    o = diff(s(row_orig(i), :));
    orig = [orig o > 0.5];
    
    % append hat
    h = diff(s_hat(row_hat(i), :));
    r = corr(s(row_orig(i), :)', s_hat(row_hat(i), :)');
    if r < 0
        h = h .* -1;
    end
    hat = [hat h];
end

% plot
%figure;
%plotroc(orig, hat);

% roc
[tpr, fpr, th] = roc(orig, hat);
auc = trapz(fpr, tpr);

end
