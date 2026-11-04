# QMeku comparison with IC50

## Outcome

Exact PDB/ligand IC50 annotations were available for 3 of 12 complexes. This is too small for validation or inferential correlation; the coefficients below are descriptive diagnostics only.

Replicate IC50 values were combined with a geometric mean because IC50 is analyzed on a logarithmic concentration scale. All individual values remain in `ic50_measurements.csv`.

| PDB | Ligand | IC50 values (nM) | Geometric mean (nM) | pIC50 | xTB | MOPAC | Sequence identity |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1ADL | ACD | 2750, 5130 | 3.76e+03 | 5.425 | -6.486 | -0.724 | 91% |
| 1BCU | PRL | 12000 | 1.2e+04 | 4.921 | -21.941 | -11.950 | 100% |
| 5F63 | 5W2 | 9.55, 1.3 | 3.52 | 8.453 | -40.934 | -30.168 | 99% |

## Descriptive correlations

- xTB/ALPB: Pearson r=-0.827; Spearman ρ=-0.500; n=3.
- MOPAC/COSMO: Pearson r=-0.868; Spearman ρ=-0.500; n=3.

A favorable association is negative because stronger potency raises pIC50 while stronger computed interaction makes the energy more negative. With only three heterogeneous targets, neither the coefficient nor its p-value is suitable as publication evidence.

## Complexes without exact IC50 annotation

2V00/V15, 3P5O/EAM, 3COW/52H, 1AMW/ADP, 1BM7/FLF, 4WIV/3P2, 1E66/HUX, 2HB1/512, 3PTB/BEN.

## Interpretation constraints

- IC50 is assay-dependent and is not an equilibrium binding constant. pIC50 is the defensible comparison scale.
- `ic50_dg_equivalent_kcal_mol` is supplied only as a plotting convenience (`RT ln(IC50 in M)`); it is not an experimental binding free energy.
- Converting IC50 to Ki requires the assay mechanism and quantities such as substrate/tracer concentration and Km or Kd.
- 1ADL has only 91% reference-sequence identity in the BindingDB annotation and should be sensitivity-flagged.
- The three observations span unrelated targets and special structural contexts; no global calibration should be fitted.

## Recommendation for a publication

Present this IC50 analysis as a secondary, explicitly exploratory panel. For a defensible IC50 benchmark, expand the crystal set prospectively using one assay family per target, exact ligand identity, documented construct/sequence, and preferably at least 20–30 independent complexes. Keep Kd, Ki, and IC50 in separate strata rather than converting and pooling them silently.

## Provenance

IC50 values were read from each downloaded RCSB `entry.json` field `rcsb_binding_affinity`, whose provenance is BindingDB. Raw rows, including source URLs and sequence identity, are in `ic50_measurements.csv`.
