function [auc, pd, pf] = roc_auc(labels, scores)
%ROC_AUC Compute ROC curve and area under curve.

labels = labels(:) > 0;
scores = scores(:);
valid = ~isnan(scores);
labels = labels(valid);
scores = scores(valid);

[scores, order] = sort(scores, 'descend');
labels = labels(order);

positive = sum(labels);
negative = numel(labels) - positive;
if positive == 0 || negative == 0
    error('Ground-truth labels must contain both anomaly and background pixels.');
end

tp = cumsum(labels);
fp_count = cumsum(~labels);

change = [true; diff(scores) ~= 0];
pd = [0; tp(change) / positive; 1];
pf = [0; fp_count(change) / negative; 1];

[pf, unique_idx] = unique(pf, 'stable');
pd = pd(unique_idx);
auc = trapz(pf, pd);
end
