# Neutral BRD4-BD1 Crystal Benchmark Series with Experimental IC50

This directory contains 17 co-crystal structures of the **first bromodomain of human BRD4 (BRD4-BD1)** bound to **neutral small-molecule inhibitors ($q = 0$)** with documented experimental $\text{IC}_{50}$ values.

None of these structures overlap with the 3 BRD4 complexes in `benchmarks/affinity_crystal/` (`3P5O`, `5F63`, `4WIV`).

## Organization

The benchmark is organized into three strata:

1. **`pfi1_series`** (4 complexes, Fish, Filippakopoulos et al., *J. Med. Chem.* 2012, DOI: [10.1021/jm3010515](https://doi.org/10.1021/jm3010515)):
   - Congeneric optimization series from fragment hit to chemical probe PFI-1.
   - Homogeneous assay format: AlphaScreen, BRD4-BD1.
   - Potency range: $220\text{ nM}$ to $23\ \mu\text{M}$ ($42\times$ dynamic range).
   - Cases: `4E96` (`0NS`), `4HBY` (`13F`), `4HBX` (`14X`), `4HBV` (`15E`).

2. **`wpf_shelf_series`** (3 complexes, Lucas, Bamborough et al., *Bioorg. Med. Chem.* 2018, DOI: [10.1016/j.bmc.2018.05.003](https://doi.org/10.1016/j.bmc.2018.05.003)):
   - Congeneric dimethylisoxazole series probing the WPF shelf region.
   - Ultra-high resolution: $1.28$ to $1.34\text{ \AA}$.
   - Homogeneous assay format: AlphaScreen, BRD4-BD1.
   - Cases: `6FSY` (`E5Q`), `6FT4` (`E5W`), `6FT3` (`E5T`).

3. **`high_res_neutral`** (10 complexes):
   - Diverse neutral chemotypes across sub-angstrom ($0.95\text{ \AA}$) to high resolution ($\le 1.63\text{ \AA}$).
   - Covers wide potency span from $251\text{ nM}$ ($\text{pIC}_{50} = 6.60$) to $126\ \mu\text{M}$ ($\text{pIC}_{50} = 3.90$).
   - Cases: `6FO5` ($0.95\text{ \AA}$), `7MRA` ($1.16\text{ \AA}$), `7MR9` ($1.19\text{ \AA}$), `4NUD` ($1.20\text{ \AA}$), `4PCE` ($1.29\text{ \AA}$), `4NUC` ($1.40\text{ \AA}$), `4NUE` ($1.30\text{ \AA}$), `4A9L` ($1.60\text{ \AA}$), `4UYD` ($1.37\text{ \AA}$), `5MLI` ($1.63\text{ \AA}$).

## Files in each case directory

```text
<set>/<PDBID>/
├── raw/
│   ├── structure.pdb
│   ├── structure.cif
│   ├── entry.json
│   ├── <CCD>_ideal.sdf
│   ├── ligand_crystal.pdb
│   ├── ligand_crystal_unprepared.sdf
│   └── receptor_protein_only_unprepared.pdb
└── metadata.json
```

All 17 cases have verified non-zero raw coordinates and exact RDKit-matched ligand SDFs preserving the crystallographic coordinates.

## QMEku benchmark results

All 17 complexes were prepared and calculated with QMeku 1.0 using both xTB and MOPAC. The reproducible workflow is:

```bash
python benchmarks/brd4_neutral_series/prepare_benchmark.py
python benchmarks/brd4_neutral_series/run_benchmark.py --method both --workers 2
python benchmarks/brd4_neutral_series/collect_results.py
python benchmarks/brd4_neutral_series/analyze_results.py
```

The complete energies are in [`results.csv`](results.csv), and the statistical interpretation is documented in [`analysis/report.md`](analysis/report.md). The primary global result is a Spearman correlation of `rho = -0.723` between MOPAC binding energy and experimental pIC50 (`n = 17`, permutation `p = 0.00143`, Holm-adjusted `p = 0.0114`).
