function alpha = solve_weight(V, index, neighbors, num_neighbors)
%SOLVE_WEIGHT Collaborative representation weights for one superpixel.

superpixel_data = V(3:end, :);
reduced = pca_reduce(superpixel_data, num_neighbors - 1);
neighbor_data = reduced(:, neighbors);
target = reduced(:, index);

G = zeros(num_neighbors - 1, num_neighbors);
h = zeros(num_neighbors, 1);
for i = 1:num_neighbors
    G(:, i) = compute_bi(neighbor_data, i, num_neighbors);
    h(i) = max(G(:, i)' * neighbor_data);
end

if num_neighbors == 2
    denom = h - (G .* neighbor_data)';
else
    denom = h - sum(G .* neighbor_data)';
end
alpha = (h - G' * target) ./ max(abs(denom), eps) .* sign(denom);

if all(alpha >= 0)
    alpha = alpha / max(sum(alpha), eps);
    return;
end

x = mean(neighbor_data, 2);
d0 = zeros(num_neighbors - 1, 1);
Dj = zeros(num_neighbors - 1, num_neighbors);
Vj = zeros(num_neighbors - 1, num_neighbors);
mu = 1 / max(norm(neighbor_data), eps);
max_iter = 100;

for iter = 1:max_iter
    x_old = x;
    v0 = (target + mu * (x - d0)) / (mu + 1);
    for j = 1:num_neighbors
        nj = x - Dj(:, j);
        if G(:, j)' * nj <= h(j)
            Vj(:, j) = nj;
        else
            Vj(:, j) = nj - ((G(:, j)' * nj - h(j)) / max(norm(G(:, j))^2, eps)) * G(:, j);
        end
    end
    x = (v0 + d0 + sum(Vj + Dj, 2)) / (num_neighbors + 1);
    d0 = d0 - x + v0;
    Dj = Dj - repmat(x, 1, num_neighbors) + Vj;

    pri_res = sqrt(norm(x - v0, 'fro')^2 + norm(repmat(x, 1, num_neighbors) - Vj, 'fro')^2);
    pri_tol = 1e-3 * max( ...
        sqrt(norm(v0, 'fro')^2 + norm(Vj, 'fro')^2), ...
        sqrt(norm(x, 'fro')^2 + norm(repmat(x, 1, num_neighbors), 'fro')^2));
    dual_res = mu * sqrt(norm(repmat(x - x_old, 1, num_neighbors + 1), 'fro')^2);
    dual_tol = 1e-3 * mu * norm([d0, Dj], 'fro');
    if pri_res <= pri_tol && dual_res <= dual_tol
        break;
    end
end

if num_neighbors == 2
    denom = h - (G .* neighbor_data)';
else
    denom = h - sum(G .* neighbor_data)';
end
alpha = (h - G' * x) ./ max(abs(denom), eps) .* sign(denom);
alpha(alpha < 0) = 0;
alpha = alpha / max(sum(alpha), eps);
end

function bi = compute_bi(points, i, num_points)
other = setdiff(1:num_points, i);
other_points = points(:, other);
base = other_points(:, end);
if num_points > 2
    A = other_points(:, 1:end-1) - base * ones(1, num_points - 2);
    bi = (eye(num_points - 1) - A * pinv(A' * A) * A') * (base - points(:, i));
else
    bi = base - points(:, i);
end
bi = bi / max(norm(bi), eps);
end
