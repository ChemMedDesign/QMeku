# Cases with a positive MOPAC energy

This folder gathers copies of every complex found under `benchmarks/` whose `mopac_e_bind` result is greater than zero. The original files stay where they are and were not modified.

Each subfolder contains, where available:

- `prepared/receptor.pdb`: the prepared protein used by QMeku;
- `prepared/ligand.sdf`: the corresponding ligand;
- `prepared/preparation.json` and `pdb4amber.log`: the preparation trail;
- `results/poses.csv`: the original result;
- `metadata.json`: metadata and provenance of the crystal.

The `raw/` reference crystal files are not distributed with this repository; they are recovered with the `download_crystals.py` script of the set of origin. The `results/qmeku.log` of each run is not distributed either: every value it prints is a column of `poses.csv`. `manifest.csv` still lists both among the files that `collect_mopac_positive.py` copies when they are present.

## Cases gathered

| PDB | xTB (kcal/mol) | MOPAC (kcal/mol) | Set of origin |
|---|---:|---:|---|
| 3COW | -10.524 | +46.730 | affinity_crystal/core_absolute |
| 1AMW | -32.729 | +23.771 | affinity_crystal/core_absolute |
| 4NUD | -2.369 | +23.186 | brd4_neutral_series/high_res_neutral |
| 2V00 | -19.357 | +17.187 | affinity_crystal/core_absolute |
| 3MXF | -26.204 | +12.622 | brd4_expanded/kd_itc_series |

`manifest.csv` records the original path, the charges, the size of the pocket and the energies. To update the selection after new calculations:

```bash
python benchmarks/collect_mopac_positive.py
```

A positive energy does not by itself imply an execution error: it means that, for that preparation, `E(complex) - E(pocket) - E(ligand)` was unfavourable with PM6-D3H4X/COSMO. The desolvation, the protonation and tautomerism, the crystallographic waters, the contacts and the sensitivity to the size of the pocket all need to be examined.

The mechanistic analysis and the comparison with controls are in [`analysis/report.md`](analysis/report.md). The structured data are recompiled with:

```bash
python benchmarks/mopac_positive_cases/analyze_cases.py
```
