# Dataset Placement

Place datasets here as MATLAB `.mat` files.

The demo expects a file with:

- `data`: hyperspectral cube, size `H x W x B`;
- `map`: binary ground-truth anomaly mask, size `H x W` (optional, but needed for AUC reporting).

Example:

```matlab
load(fullfile('data', 'San_Diego.mat'))  % should load data and map
```

Datasets are not bundled in this public folder to avoid redistributing data with unclear or restricted licenses.
