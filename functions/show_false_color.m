function show_false_color(data)
%SHOW_FALSE_COLOR Display a simple false-color view of a hyperspectral cube.

[~, ~, bands] = size(data);
if bands >= 30
    idx = unique(max(1, min(bands, [30, 20, 10])), 'stable');
elseif bands >= 3
    idx = [bands, ceil(bands / 2), 1];
else
    idx = ones(1, 3);
end

rgb = data(:, :, idx);
if size(rgb, 3) == 1
    rgb = repmat(rgb, 1, 1, 3);
elseif size(rgb, 3) == 2
    rgb(:, :, 3) = rgb(:, :, 2);
end
imshow(normalize_minmax(rgb));
end
