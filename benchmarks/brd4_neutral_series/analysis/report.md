# QMeku benchmark on neutral BRD4-BD1 crystals

Date of the analysis: 2026-10-31 (America/Santiago)

## Executive summary

**17/17 crystal complexes** of BRD4-BD1 with neutral ligands were prepared, validated and computed. Every model used by QMeku had a ligand charge of 0, an even number of electrons and complete results with both xTB and MOPAC. There were no calculation failures.

The clearest statistical result was the global monotonic correlation of MOPAC with the experimental activity: **Spearman ρ = −0.723**, bootstrap 95% CI [−0.876, −0.382], permutation *p* = 0.00143 and Holm-adjusted *p* = 0.0114. The negative sign is the expected one: a more negative computed energy corresponds to a higher pIC50. For xTB the global association was smaller and inconclusive: **ρ = −0.397**, 95% CI [−0.752, 0.071], *p* = 0.116.

These results support a **ranking ability** for MOPAC on this set, but not a calibration of absolute energy nor a reproduction of IC50. The global inference is exploratory because it pools several assays and chemotypes.

## Preparation and quality control

- Input: the co-crystallised ligand and the BRD4-BD1 receptor of each structure.
- Receptor: prepared with `pdb4amber --reduce`, hydrogens added and the crystallographic heavy-atom coordinates kept.
- Ligand: neutral protomer derived from the CCD, crystallographic heavy-atom coordinates preserved and hydrogens added with RDKit.
- Calculation: QMeku 1.0, implicit water, `--total-only --no-gas`, xTB and MOPAC, without the per-residue decomposition.
- Validation: 17 ligands with formal charge 0; 17 models with an even number of electrons; 17 complete result files.

## Computed energies

| PDB | Series | pIC50 | xTB (kcal/mol) | MOPAC (kcal/mol) |
|---|---|---:|---:|---:|
| 4E96 | PFI-1 | 6.261 | −32.067 | −17.297 |
| 4HBY | PFI-1 | 5.357 | −27.397 | −13.440 |
| 4HBX | PFI-1 | 5.319 | −28.834 | −15.966 |
| 4HBV | PFI-1 | 4.638 | −15.070 | −11.793 |
| 6FSY | WPF shelf | 6.780 | −20.011 | −20.725 |
| 6FT4 | WPF shelf | 6.630 | −20.789 | −20.710 |
| 6FT3 | WPF shelf | 6.422 | −17.405 | −17.888 |
| 6FO5 | high resolution | 6.600 | −38.488 | −25.390 |
| 7MRA | high resolution | 6.548 | −36.865 | −21.238 |
| 7MR9 | high resolution | 6.000 | −31.575 | −24.688 |
| 4NUD | high resolution | 5.365 | −2.369 | +23.186 |
| 4PCE | high resolution | 5.155 | −22.274 | −18.221 |
| 4NUC | high resolution | 4.565 | −15.340 | −4.428 |
| 4NUE | high resolution | 4.527 | −21.337 | −13.352 |
| 4A9L | high resolution | 4.523 | −25.992 | −13.726 |
| 4UYD | high resolution | 4.100 | −15.849 | −10.512 |
| 5MLI | high resolution | 3.900 | −16.210 | −12.670 |

pIC50 was recomputed without rounding as `9 − log10(IC50[nM])`, using the geometric mean documented in the manifest when there were several measurements.

## Statistical analysis

Spearman was the primary analysis because the aim is to test monotonic ranking, not identity between scales. Pearson was secondary and assesses linearity. The Spearman *p* values are exact for `n ≤ 8` and were obtained with 100,000 permutations for larger panels. The 95% CIs are percentiles of 20,000 paired resamples. The Holm correction covered the eight Spearman comparisons (four panels × two methods).

### Primary monotonic correlation

| Panel | n | xTB ρ (95% CI) | *p* | MOPAC ρ (95% CI) | *p* |
|---|---:|---:|---:|---:|---:|
| PFI-1, homogeneous assay | 4 | −0.800 [−1.000, 1.000] | 0.333 | −0.800 [−1.000, 1.000] | 0.333 |
| WPF shelf, homogeneous assay | 3 | −0.500 [−1.000, 1.000] | 1.000 | −1.000 [−1.000, −1.000] | 0.333 |
| High resolution, heterogeneous | 10 | −0.564 [−0.963, 0.250] | 0.096 | −0.588 [−0.962, 0.245] | 0.081 |
| All | 17 | −0.397 [−0.752, 0.071] | 0.116 | **−0.723 [−0.876, −0.382]** | **0.00143** |

After Holm, only MOPAC on the global panel stayed below 0.05 (adjusted *p* = 0.0114). The intervals of the panels with `n = 3–4` are deliberately wide: a large correlation there is not by itself stable evidence.

### Pearson and sensitivity to influential points

| Panel | xTB r (*p*) | MOPAC r (*p*) |
|---|---:|---:|
| PFI-1 | −0.892 (0.108) | −0.890 (0.110) |
| WPF shelf | −0.795 (0.415) | −0.911 (0.271) |
| High resolution | −0.661 (0.0375) | −0.315 (0.376) |
| All | −0.464 (0.0606) | −0.409 (0.103) |

The nominal Pearson result for xTB on the high-resolution panel is not robust: its bootstrap 95% CI is [−0.984, 0.107]. It must therefore not be presented as a confirmed linear relation.

**4NUD is the most influential point.** Its energy is −2.369 kcal/mol with xTB and +23.186 kcal/mol with MOPAC, despite a pIC50 of 5.365. Omitting it as a diagnostic only:

- high resolution: xTB `r = −0.939`, MOPAC `r = −0.841`;
- complete set: xTB `r = −0.551`, MOPAC `r = −0.801` and MOPAC `ρ = −0.806`.

4NUD was not removed from the primary analysis: there was no execution error, no incorrect charge and no odd-electron model. Its effect explains why MOPAC shows a clear rank correlation but a more modest global linear one. Its protonation state, water network and pocket contacts should be examined specifically before its absolute energy is interpreted.

### Molecular size

No conclusive global relation was found between the energy and the number of heavy atoms (xTB `ρ = −0.358`, *p* = 0.159; MOPAC `ρ = −0.210`, *p* = 0.419). Dividing the energy by the number of heavy atoms, the global association with pIC50 was small (xTB `ρ = 0.115`; MOPAC `ρ = −0.174`). This normalisation is diagnostic only and is no substitute for a physically validated energy function.

## Interpretation and limits

1. **MOPAC ranks this dataset better than xTB.** The advantage is seen on the global panel and is directionally consistent in PFI-1, WPF shelf and high resolution.
2. **The evidence within each series is still limited by sample size.** PFI-1 (`n = 4`) and WPF (`n = 3`) show useful trends, but they do not allow precise estimates.
3. **IC50 is not a binding free energy.** It depends on the assay format, the concentrations, the construct and the experimental conditions. It must not be turned into an expectation of numerical equality with QMeku.
4. **The global panel mixes assays and chemotypes.** Its correlation may include differences between series; it should be treated as exploratory until the series with homogeneous measurements are enlarged.
5. **A crystal structure is a single conformation.** Conformational entropy, dynamics, conserved explicit waters and multiple microstates were not included.

The appropriate conclusion for publication is that, on this neutral panel, **QMeku/MOPAC shows a monotonic association with the experimental activity**, whereas **QMeku/xTB gives no conclusive global evidence**. Neither quantitative prediction of IC50 nor causality is demonstrated.

## Reproducibility

```bash
python benchmarks/brd4_neutral_series/prepare_benchmark.py
python benchmarks/brd4_neutral_series/run_benchmark.py --method both --workers 2
python benchmarks/brd4_neutral_series/collect_results.py
python benchmarks/brd4_neutral_series/analyze_results.py
```

Files produced:

- `results.csv`: activity, control descriptors and energies per complex.
- `analysis/correlation_metrics.csv`: correlations, intervals, tests and leave-one-out sensitivity.
- `analysis/descriptives.csv`: ranges per panel.
- `analysis/size_diagnostics.csv`: molecular-size diagnostics.
- `analysis/correlation_panels.png` and `.pdf`: correlation figure.

The statistical design followed the recommendations of **Scientific Agent Skills**: Kassis T. et al., *Scientific Agent Skills: A Skill Library for Scientific Agents*, arXiv:2609.00065v2, DOI: [10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065).
