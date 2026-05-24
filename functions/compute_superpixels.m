function [labels, num_segments] = compute_superpixels(img, superpixel_size, compactness)
%COMPUTE_SUPERPIXELS SLIC wrapper with practical fallbacks.

if nargin < 3 || isempty(compactness)
    compactness = 0.3;
end

if exist('vl_slic', 'file') == 3 || exist('vl_slic', 'file') == 2
    norm_img = sqrt(sum(img.^2, 3));
    scale = mean(norm_img(:));
    if scale > 0
        img = img / scale;
    end
    labels = double(vl_slic(single(img), superpixel_size, compactness)) + 1;
    num_segments = max(labels(:));
    return;
end

if exist('superpixels', 'file') == 2
    rgb = normalize_minmax(img);
    if size(rgb, 3) == 1
        rgb = repmat(rgb, 1, 1, 3);
    elseif size(rgb, 3) > 3
        rgb = rgb(:, :, 1:3);
    end
    target_segments = max(2, round(numel(rgb(:, :, 1)) / max(superpixel_size^2, 1)));
    labels = double(superpixels(rgb, target_segments, 'Compactness', max(1, compactness * 20)));
    num_segments = max(labels(:));
    return;
end

warning('CWP_RPCA:GridSuperpixels', ...
    ['Neither vl_slic nor MATLAB superpixels was found. ' ...
    'Using a grid fallback; install VLFeat for paper-faithful runs.']);
[height, width, ~] = size(img);
labels = zeros(height, width);
label = 0;
for r = 1:superpixel_size:height
    for c = 1:superpixel_size:width
        label = label + 1;
        rr = r:min(r + superpixel_size - 1, height);
        cc = c:min(c + superpixel_size - 1, width);
        labels(rr, cc) = label;
    end
end
num_segments = label;
end
