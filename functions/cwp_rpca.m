function [score_map, elapsed, history, model] = cwp_rpca(data, B_csrd, params)
%CWP_RPCA Confidence-weighted prior-guided RPCA.
%
% Inputs
%   data   : H x W x B normalized hyperspectral cube.
%   B_csrd : B x (H*W) CSRD background prior.
%   params : structure from default_cwp_params().
%
% Outputs
%   score_map : H x W anomaly score map in [0, 1].
%   elapsed   : optimization time in seconds.
%   history   : ADMM objective and residual history.
%   model     : selected internal variables useful for diagnostics.

[height, width, bands] = size(data);
num_pixels = height * width;
Y = reshape(data, num_pixels, bands)';

lambda_a = get_param(params, 'lambda_a', 0.05);
lambda_b = get_param(params, 'lambda_b', 1e-3);
lambda_tv = get_param(params, 'lambda_tv', 1e-3);
rank_n = get_param(params, 'N', 1);
threshold = get_param(params, 'threshold', 0.8);
alpha = get_param(params, 'alpha', 2);
max_iter = get_param(params, 'max_iter', 100);
rho = get_param(params, 'rho', 1.1);
tol = get_param(params, 'tol', 1e-3);

S = zeros(rank_n, num_pixels);
A = zeros(bands, num_pixels);

U1 = zeros(rank_n, num_pixels); J1 = zeros(rank_n, num_pixels);
U2 = zeros(bands, num_pixels);  J2 = zeros(bands, num_pixels);
U3 = zeros(bands, num_pixels);  J3 = zeros(bands, num_pixels);
U4 = zeros(rank_n, num_pixels); J4 = zeros(rank_n, num_pixels);
U5 = zeros(rank_n, num_pixels, 2); J5 = zeros(rank_n, num_pixels, 2);

mu = 1 / max(norm(Y, 'fro'), eps);
mu_bar = mu * 1e7;

history.obj = [];
history.pri_res = [];
history.dual_res = [];
history.pri_tol = [];
history.dual_tol = [];

timer = tic;

E = compute_basis(B_csrd, rank_n);
S_prior = E' * B_csrd;

residual = sqrt(sum((Y - B_csrd).^2, 1));
residual = residual - min(residual);
if max(residual) > 0
    residual = residual / max(residual);
end
confidence = exp(-alpha * residual.^2);
lambda_b_map = lambda_b * confidence;
lambda_b_mat = repmat(lambda_b_map, rank_n, 1);

otf_dx = local_psf2otf([1, -1], [height, width]);
otf_dy = local_psf2otf([1; -1], [height, width]);
denom_kernel = 2 + abs(otf_dx).^2 + abs(otf_dy).^2;

EI = [E, eye(bands)];
G = EI' * EI + eye(bands + rank_n);
EIY = EI' * Y;

for iter = 1:max_iter
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;
    U4_old = U4;
    U5_old = U5;

    U12 = G \ (EIY + [S - J1; A - J2]);
    U1 = U12(1:rank_n, :);
    U2 = U12(rank_n + 1:end, :);

    U3 = soft_threshold(A - J3, lambda_a / mu);
    U4 = (lambda_b_mat .* S_prior + mu * (S - J4)) ./ (lambda_b_mat + mu);

    grad_S = gradient_periodic(S, height, width);
    U5 = soft_threshold(grad_S - J5, lambda_tv / mu);

    rhs = U1 + J1 + U4 + J4 - divergence_periodic(U5 + J5, height, width);
    for n = 1:rank_n
        rhs_img = reshape(rhs(n, :), height, width);
        S_img = real(ifft2(fft2(rhs_img) ./ denom_kernel));
        S(n, :) = S_img(:)';
    end

    A = (U2 + J2 + U3 + J3) / 2;

    J1 = J1 - S + U1;
    J2 = J2 - A + U2;
    J3 = J3 - A + U3;
    J4 = J4 - S + U4;
    J5 = J5 - grad_S + U5;

    data_term = 0.5 * norm(Y - E * S - A, 'fro')^2;
    sparse_term = lambda_a * sum(abs(A(:)));
    prior_diff = S - S_prior;
    prior_term = 0.5 * sum(lambda_b_map .* sum(prior_diff.^2, 1));
    tv_term = lambda_tv * sum(abs(grad_S(:)));
    history.obj(end + 1) = data_term + sparse_term + prior_term + tv_term;

    pri_res = sqrt(norm(S - U1, 'fro')^2 + norm(A - U2, 'fro')^2 + ...
        norm(A - U3, 'fro')^2 + norm(S - U4, 'fro')^2 + ...
        norm(grad_S(:) - U5(:), 2)^2);
    pri_tol = tol * max( ...
        sqrt(norm([U1, U4], 'fro')^2 + norm(U5(:), 2)^2), ...
        sqrt(norm([S, S], 'fro')^2 + norm(grad_S(:), 2)^2));

    dual_res = mu * sqrt(norm(U1 - U1_old, 'fro')^2 + ...
        norm(U2 - U2_old, 'fro')^2 + norm(U3 - U3_old, 'fro')^2 + ...
        norm(U4 - U4_old, 'fro')^2 + norm(U5(:) - U5_old(:), 2)^2);
    dual_tol = tol * mu * sqrt(norm([J1, J4], 'fro')^2 + ...
        norm([J2, J3], 'fro')^2 + norm(J5(:), 2)^2);

    history.pri_res(end + 1) = pri_res;
    history.dual_res(end + 1) = dual_res;
    history.pri_tol(end + 1) = pri_tol;
    history.dual_tol(end + 1) = dual_tol;

    if pri_res <= pri_tol && dual_res <= dual_tol
        break;
    end

    mu = min(mu * rho, mu_bar);
end

score = sqrt(sum(A.^2, 1));
score_map = reshape(score, height, width);
score_map = normalize_minmax(score_map);
score_map = cs_transform(score_map, threshold);

elapsed = toc(timer);

model.E = E;
model.S = S;
model.A = A;
model.S_prior = S_prior;
model.confidence = reshape(confidence, height, width);
model.lambda_b_map = reshape(lambda_b_map, height, width);
end

function value = get_param(params, name, default_value)
if isfield(params, name)
    value = params.(name);
else
    value = default_value;
end
end

function y = soft_threshold(x, lambda)
y = sign(x) .* max(abs(x) - lambda, 0);
end

function G = gradient_periodic(S, height, width)
[rank_n, num_pixels] = size(S);
G = zeros(rank_n, num_pixels, 2);
for n = 1:rank_n
    img = reshape(S(n, :), height, width);
    gx = circshift(img, [0, -1]) - img;
    gy = circshift(img, [-1, 0]) - img;
    G(n, :, 1) = gx(:)';
    G(n, :, 2) = gy(:)';
end
end

function D = divergence_periodic(V, height, width)
[rank_n, num_pixels, ~] = size(V);
D = zeros(rank_n, num_pixels);
for n = 1:rank_n
    vx = reshape(V(n, :, 1), height, width);
    vy = reshape(V(n, :, 2), height, width);
    dx = vx - circshift(vx, [0, 1]);
    dy = vy - circshift(vy, [1, 0]);
    D(n, :) = reshape(dx + dy, 1, []);
end
end

function otf = local_psf2otf(psf, out_size)
pad = zeros(out_size);
psf_size = size(psf);
pad(1:psf_size(1), 1:psf_size(2)) = psf;
for dim = 1:numel(psf_size)
    pad = circshift(pad, -floor(psf_size(dim) / 2), dim);
end
otf = fft2(pad);
end
