function E = compute_basis(X, rank_n)
%COMPUTE_BASIS Leading left singular vectors of a background prior.

[U, ~, ~] = svd(X, 'econ');
rank_n = min(rank_n, size(U, 2));
E = U(:, 1:rank_n);
end
