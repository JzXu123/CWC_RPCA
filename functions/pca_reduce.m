function [scores, coeff] = pca_reduce(X, n_components)
%PCA_REDUCE PCA scores for columns of X.

[bands, num_samples] = size(X);
n_components = min([n_components, bands, num_samples]);
mu = mean(X, 2);
X0 = X - mu;
[U, ~, ~] = svd(X0, 'econ');
coeff = U(:, 1:n_components);
scores = coeff' * X0;
end
