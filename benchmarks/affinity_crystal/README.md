# Crystal-affinity benchmark

This directory records crystal complexes selected for a prospective comparison of
QMeku total interaction energies with experimental affinity labels.

## Fixed calculation protocol

- QMeku 1.0
- crystallographic ligand pose; no geometry optimisation
- 4.0 A per-pose pocket
- GFN2-xTB/ALPB water and PM6-D3H4X/COSMO water
- `--total-only --no-gas`
- no per-residue decomposition

The calculated `E_bind` is a single-point ligand-pocket interaction energy, not a
standard binding free energy. Experimental `Kd` labels are converted to standard
free energy at 298.15 K only for comparison. `Ki` cases are flagged and should be
analysed separately or treated as approximate equivalents.

## State of the set

All 12 complexes have validated `prepared/` inputs and results:

- **All 12 complexes were computed with QMeku 1.0** using the fixed protocol above.
  Their `poses.csv`, the consolidated `results.csv`, the correlation analysis and the
  IC50 comparison all come from that one run. The readable `qmeku.log` of each run is
  not distributed: every value it prints is a column of `poses.csv`.
- `1BCU` (thrombin) was re-prepared because its PDB insertion codes (`36A`, `184A`,
  `186A`–`186D`, `221A`) previously caused `pdb4amber --reduce --add-missing-atoms`
  to create spurious internal termini. The corrected receptor changes the pocket from
  253 atoms/charge -6 to 248 atoms/charge -5 and changes its energies from xTB -61.550
  and MOPAC -8.376 to xTB -21.941 and MOPAC -11.950 kcal/mol. The results stored here
  are those of the corrected receptor.
- The `raw/` directories of the downloaded crystals and the wwPDB validation reports
  are **not distributed** with this repository: they are reproducible with
  `download_crystals.py`, and `metadata.json` records the SHA-256 of every downloaded
  file, so a new download can be checked against the original one.

## Download the deposited crystals

The downloader performs a bounded retrieval of the 11 pending PDB entries from the
official RCSB endpoints:

```bash
python benchmarks/affinity_crystal/download_crystals.py --dry-run
python benchmarks/affinity_crystal/download_crystals.py
```

It downloads the deposited PDB and mmCIF files, RCSB entry metadata, the CCD ideal
ligand SDF, and the wwPDB validation report. For 1BM7 it additionally downloads
biological assembly 1. Each case is organized as:

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
├── prepared/
├── validation/wwpdb_validation.pdf
├── results/
└── metadata.json
```

`raw/ligand_crystal_unprepared.sdf` is generated preferentially with RDKit by
transferring bond orders and formal charges from the CCD ideal template onto the
crystal coordinates. Open Babel is used only as a fallback. Its stereochemistry,
formal charge, and protonation must still be checked before it is copied to
`prepared/ligand.sdf`.
`raw/receptor_protein_only_unprepared.pdb` has no waters, ions, or ligands, but it
is deliberately **not** called `receptor.pdb`: hydrogens, histidine states,
disulfides, termini, missing atoms, modified residues, and AMBER names still need
to be validated before creating `prepared/receptor.pdb`.

## Prepare QMeku inputs

With AmberTools and RDKit available, prepare and validate all downloaded cases:

```bash
python benchmarks/affinity_crystal/prepare_benchmark.py
```

The script preserves crystallographic heavy-atom coordinates, uses
`pdb4amber --reduce` for protein hydrogens and AMBER residue names, restores the
two modified polymer residues in 1ADL to their standard parents, and assigns the
documented pH 7.4 ligand protomers. It writes `prepared/preparation.json` and
constructs the same 4 A pocket that QMeku will use as an input validation.

Geometric chain breaks are marked before protonation so AmberTools creates valid
termini instead of leaving odd-electron fragments. Missing protein side-chain atoms
are rebuilt with AmberTools; this is material for 1E66 and is recorded in its
`pdb4amber.log`. Existing ligand heavy atoms are never optimized: every prepared
ligand records a 0.0 A heavy-atom RMSD from the downloaded crystal pose.

Run all cases whose prepared inputs exist:

```bash
python benchmarks/affinity_crystal/run_benchmark.py
python benchmarks/affinity_crystal/collect_results.py
```

## Correlation analysis

After all energies are collected, run:

```bash
python benchmarks/affinity_crystal/analyze_correlations.py
```

The analysis keeps `Kd`, `Ki`, and the three-ligand BRD4 ranking separate. It
reports Pearson and Spearman coefficients, paired 95% bootstrap intervals for the
`Kd` panels (20,000 resamples; seed 20261031), leave-one-out sensitivity, and
publication-oriented PNG/PDF figures under `analysis/`. No observations are
excluded after inspecting the results.

Prepared inputs must be placed at the paths recorded in `manifest.csv`. Receptors
must contain hydrogens, AMBER residue names, and amino acids only. Ligand SDF files
must preserve the exact crystallographic coordinates, bond orders, stereochemistry,
formal charge, and protonation decision.

## IC50 comparison

The downloaded RCSB `entry.json` files include BindingDB affinity annotations when
available. Extract exact PDB/ligand IC50 measurements, aggregate replicate values
in log space, and generate the exploratory comparison with:

```bash
python benchmarks/affinity_crystal/analyze_ic50.py
```

Outputs are written under `analysis_ic50/`. IC50 is compared as pIC50; the script
does not treat IC50 as a thermodynamic binding constant or silently pool it with
Kd/Ki. The current exact-annotation coverage is only 3/12 complexes, so these
statistics are descriptive rather than a validation result.
