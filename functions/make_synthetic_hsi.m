function [data, gt_map] = make_synthetic_hsi()
%MAKE_SYNTHETIC_HSI Small synthetic HSI for smoke testing the public demo.

height = 40;
width = 40;
bands = 16;
[x_grid, y_grid] = meshgrid(linspace(0, 1, width), linspace(0, 1, height));

spectral_1 = sin(linspace(0, pi, bands));
spectral_2 = cos(linspace(0, 2 * pi, bands));
background = zeros(height, width, bands);
for b = 1:bands
    background(:, :, b) = 0.55 * spectral_1(b) * x_grid + ...
        0.35 * spectral_2(b) * y_grid + 0.08 * randn(height, width);
end

gt_map = false(height, width);
gt_map(9:12, 27:30) = true;
gt_map(28:31, 12:15) = true;

anomaly_signature = reshape(linspace(0.2, 1.2, bands), 1, 1, bands);
data = background;
for b = 1:bands
    band = data(:, :, b);
    band(gt_map) = band(gt_map) + anomaly_signature(:, :, b);
    data(:, :, b) = band;
end

data = normalize_minmax(data);
end
