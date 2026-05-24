function params = default_cwp_params()
%DEFAULT_CWP_PARAMS Default parameters for the public CWP-RPCA demo.

params.lambda_a = 0.05;
params.lambda_b =1e-4;
params.lambda_tv = 1e-1;
params.N = 10;
params.threshold = 0.8;
params.alpha = 2;

params.superpixel_size = 15;
params.compactness = 0.3;

params.max_iter = 100;
params.rho = 1.1;
params.tol = 1e-3;
end
