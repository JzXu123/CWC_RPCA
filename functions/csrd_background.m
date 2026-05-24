function [B_csrd, elapsed, labels] = csrd_background(data, superpixel_size, compactness)
%CSRD_BACKGROUND Collaborative superpixel representation background prior.

timer = tic;
[height, width, bands] = size(data);
num_pixels = height * width;
Y = reshape(data, num_pixels, bands)';

pca_img = pca_reduce(Y, min(3, bands));
pca_img = reshape(pca_img', height, width, []);
[labels, num_segments] = compute_superpixels(pca_img, superpixel_size, compactness);
labels = relabel_contiguous(labels);
num_segments = max(labels(:));

coords = cell(num_segments, 1);
centers_xy = zeros(2, num_segments);
centers_spec = zeros(bands, num_segments);

for k = 1:num_segments
    idx = find(labels == k);
    coords{k} = idx;
    [row, col] = ind2sub([height, width], idx);
    centers_xy(:, k) = [mean(row); mean(col)];
    centers_spec(:, k) = mean(Y(:, idx), 2);
end

V = [centers_xy; centers_spec];
neighbors = segment_neighbors(labels, num_segments);
recon_centers = zeros(bands, num_segments);

for k = 1:num_segments
    nb = neighbors{k};
    if isempty(nb)
        recon_centers(:, k) = centers_spec(:, k);
    elseif numel(nb) == 1
        recon_centers(:, k) = centers_spec(:, nb);
    else
        weights = solve_weight(V, k, nb, numel(nb));
        recon_centers(:, k) = centers_spec(:, nb) * weights;
    end
end

B_csrd = zeros(bands, num_pixels);
for k = 1:num_segments
    B_csrd(:, coords{k}) = repmat(recon_centers(:, k), 1, numel(coords{k}));
end

elapsed = toc(timer);
end

function labels = relabel_contiguous(labels)
ids = unique(labels(:));
new_labels = zeros(size(labels));
for i = 1:numel(ids)
    new_labels(labels == ids(i)) = i;
end
labels = new_labels;
end

function neighbors = segment_neighbors(labels, num_segments)
neighbors = cell(num_segments, 1);

top_labels = labels(1:end-1, :);
bottom_labels = labels(2:end, :);
left_labels = labels(:, 1:end-1);
right_labels = labels(:, 2:end);
pairs = [top_labels(:), bottom_labels(:); left_labels(:), right_labels(:)];
pairs = pairs(pairs(:, 1) ~= pairs(:, 2), :);

for i = 1:size(pairs, 1)
    a = pairs(i, 1);
    b = pairs(i, 2);
    neighbors{a}(end + 1) = b;
    neighbors{b}(end + 1) = a;
end

for k = 1:num_segments
    neighbors{k} = unique(neighbors{k});
end
end
