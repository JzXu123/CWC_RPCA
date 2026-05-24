clear; clc; close all;

addpath(genpath(fullfile(pwd, 'functions')));

dataset_file = fullfile(pwd, 'data', 'Salians.mat');
if isfile(dataset_file)
    [data, mask] = load_hsi_mat(dataset_file);
else
    fprintf(['Dataset not found: %s\n' ...
        'Running a small synthetic example instead. Add a .mat dataset ' ...
        'or edit dataset_file for real experiments.\n'], dataset_file);
    [data, mask] = make_synthetic_hsi();
end

data = normalize_minmax(double(data));

params = default_cwp_params();

fprintf('Running CSRD background prior...\n');
[B_csrd, t_csrd] = csrd_background(data, params.superpixel_size, params.compactness);

fprintf('Running CWP-RPCA...\n');
[score_map, t_cwp, history] = cwp_rpca(data, B_csrd, params);

fprintf('CSRD time: %.3f s\n', t_csrd);
fprintf('CWP-RPCA time: %.3f s\n', t_cwp);

if ~isempty(mask)
    [auc, pd, pf] = roc_auc(mask(:), score_map(:));
    fprintf('AUC: %.6f\n', auc);
else
    auc = NaN;
    pd = [];
    pf = [];
end

figure('Color', 'w', 'Name', 'CWP-RPCA demo');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
show_false_color(data);
title('False color');

nexttile;
if ~isempty(mask)
    imagesc(mask); axis image off; colormap(gca, gray);
    title('Ground truth');
else
    text(0.5, 0.5, 'No ground truth', 'HorizontalAlignment', 'center');
    axis off;
end

nexttile;
imagesc(score_map); axis image off; colormap(gca, hot); colorbar;
if isnan(auc)
    title('CWP-RPCA score');
else
    title(sprintf('CWP-RPCA score, AUC %.4f', auc));
end

if ~isempty(pd)
    figure('Color', 'w', 'Name', 'ROC');
    plot(pf, pd, 'LineWidth', 1.5);
    grid on; axis square;
    xlabel('False positive rate');
    ylabel('Detection probability');
    title(sprintf('ROC, AUC %.4f', auc));
end

if ~isempty(history.obj)
    figure('Color', 'w', 'Name', 'ADMM convergence');
    semilogy(history.obj, 'LineWidth', 1.2);
    grid on;
    xlabel('Iteration');
    ylabel('Objective');
    title('Objective history');
end
