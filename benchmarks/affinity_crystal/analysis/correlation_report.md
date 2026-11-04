# QMeku crystal-affinity correlation analysis

This is an exploratory validation. QMeku values are single-point ligand-pocket interaction energies, not standard binding free energies.

## Planned strata

- Primary cross-target panel: six `Kd` complexes in `core_absolute`; one target per complex.
- Expanded `Kd` panel: nine complexes, including repeated BRD4 measurements and the 1BCU stress case.
- BRD4 ranking: three ligands, interpreted descriptively because n=3.
- `Ki` panel: three complexes, interpreted descriptively and not pooled with `Kd`.
- 3PTB remains a development control.

## Results

### Primary cross-target Kd panel

- xTB/ALPB: Pearson r=-0.293, 95% bootstrap CI [-0.972, 0.880], p=0.573; Spearman ρ=-0.200, 95% bootstrap CI [-0.920, 1.000], p=0.704.
- MOPAC/COSMO: Pearson r=0.314, 95% bootstrap CI [-0.739, 1.000], p=0.545; Spearman ρ=0.543, 95% bootstrap CI [-0.800, 1.000], p=0.266.

### Expanded Kd panel

- xTB/ALPB: Pearson r=0.103, 95% bootstrap CI [-0.635, 0.694], p=0.791; Spearman ρ=0.200, 95% bootstrap CI [-0.583, 0.858], p=0.606.
- MOPAC/COSMO: Pearson r=0.287, 95% bootstrap CI [-0.331, 0.919], p=0.454; Spearman ρ=0.417, 95% bootstrap CI [-0.374, 0.900], p=0.265.

### BRD4 ranking

- Experimental order: 5F63 > 3P5O > 4WIV (strongest to weakest).
- xTB reproduces the exact order (Spearman ρ=1.0).
- MOPAC orders 5F63 > 4WIV > 3P5O (Spearman ρ=0.5).
- With n=3 these are descriptive ranking results, not reliable significance tests.

### Interpretation

Neither method shows evidence of a general cross-target correlation in this small heterogeneous set. The xTB result is promising for within-target BRD4 ranking but does not generalize to absolute affinity across targets. Charged systems and special contexts are influential: 1BCU and 2HB1 receive very favorable xTB interaction energies despite weak experimental affinity, while MOPAC gives positive interaction energies for 2V00, 3COW, and 1AMW.

## Diagnostics and limitations

- No observations were excluded after examining the results.
- Pearson and Spearman were both reported because linearity and normal residuals cannot be established credibly with n=6–9.
- Bootstrap intervals use paired resampling with 20,000 draws and seed 20261031; their width reflects the small sample.
- Observations in the expanded panel are not fully independent because BRD4 contributes three complexes.
- Experimental ΔG values derived from Ki are not treated as Kd measurements.
- The benchmark is too small for a stable fitted calibration model.

## Reproducibility

Source data: `../results.csv`. Analysis: `../analyze_correlations.py`.
