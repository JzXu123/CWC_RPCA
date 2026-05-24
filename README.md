# CWP-RPCA for Hyperspectral Anomaly Detection

This repository contains a clean MATLAB implementation of the method in:

> J. Xu, W. Jiang, L. Chen, C. Zhang, M. Wildgruber, X. Yang, and X. Ma, "Confidence-Weighted Prior-Guided RPCA for Hyperspectral Anomaly Detection," IEEE Signal Processing Letters, vol. 33, pp. 1766-1770, 2026. DOI: 10.1109/LSP.2026.3682998
        
        .

## Method

CWP-RPCA decomposes a hyperspectral image into a subspace low-rank background and sparse anomalies. It uses:

- a CSRD superpixel background prior;
- residual-driven confidence weights for pixel-wise prior guidance;
- a mild TV penalty on background coefficient maps;
- an FFT-accelerated ADMM solver.

## Folder Layout

```text
CWP_RPCA_public/
  demo_cwp_rpca.m          Example entry point
  functions/               Core algorithm and utilities
  data/README.md           Dataset format and placement notes
  third_party/README.md    Optional VLFeat dependency notes
```

## Quick Start

1. Optionally put a `.mat` dataset in `data/`. The file should contain:
   - `data`: an `H x W x B` hyperspectral cube;
   - `map`: an `H x W` binary ground-truth anomaly mask, if available.
2. For the original SLIC-based CSRD prior, install VLFeat and add it to the MATLAB path. See `third_party/README.md`.
3. Run:

```matlab
demo_cwp_rpca
```

If VLFeat is not available, the demo falls back to MATLAB `superpixels` when possible, and then to a grid segmentation fallback. The VLFeat path is recommended for reproducing the paper most closely.

When `data/San_Diego.mat` is not present, the demo runs on a small synthetic HSI so that the code can be tested immediately.

## Main API

```matlab
[B_csrd, t_csrd] = csrd_background(data, params.superpixel_size, params.compactness);
[score_map, t_cwp, history, model] = cwp_rpca(data, B_csrd, params);
```

`data` should be normalized to `[0, 1]` before calling the algorithm.

