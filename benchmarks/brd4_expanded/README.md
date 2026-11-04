# Expanded BRD4-BD1 crystal benchmark

This benchmark is isolated from `affinity_crystal/` and contains two internally
homogeneous seven-complex panels:

- `ic50_series`: one chemical series, BRD4-BD1 crystal structures, Oxford
  AlphaScreen assay B (four replicates) from DOI `10.1021/acschembio.7b00638`.
- `kd_itc_series`: direct BRD4-BD1 ITC measurements from DOI
  `10.1021/acs.jmedchem.1c01096`.

The calculation protocol matches the previous benchmark: the crystallographic
pose, a 4.0 A per-pose pocket, GFN2-xTB/ALPB water and PM6-D3H4X/COSMO water,
`--total-only --no-gas`, and no per-residue decomposition. QMeku energies are
single-point interaction energies and are not binding free energies.

## Reproducible workflow

```bash
python benchmarks/brd4_expanded/download_crystals.py --dry-run
python benchmarks/brd4_expanded/download_crystals.py --no-validation
python benchmarks/brd4_expanded/prepare_benchmark.py
python benchmarks/brd4_expanded/run_benchmark.py --method both --workers 2
python benchmarks/brd4_expanded/collect_results.py
python benchmarks/brd4_expanded/analyze_results.py
```

Or resume every stage with one idempotent command:

```bash
python benchmarks/brd4_expanded/run_all.py --method both --workers 2
```

Raw PDB/mmCIF coordinates, CCD ligand templates, RCSB entry metadata, SHA-256
hashes, preparation decisions, prepared inputs, logs, and energies are retained
under each PDB directory.

## Protonation contract

Ligand heavy-atom coordinates are never optimized. At pH 7.4 the seven IC50
ligands are represented as +1 cations because they contain a basic tertiary
amine. All seven ITC-panel ligands are also +1; for (+)-JQ1 the formal +1 state
is retained directly from the CCD entry (diazepin-10-ium), without selecting an
additional protonation site.
If the CCD template is neutral where +1 is required, the preparation script
protonates the preferred aliphatic tertiary amine and records all candidates and
the chosen atom in `prepared/preparation.json`. This deterministic choice should
still be treated as a protonation-state sensitivity variable, especially for the
disubstituted piperazine in volasertib.

The 6CJ1 entry contains two ligand stereoisomers. The manifest deliberately uses
F2S, which RCSB designates as the ligand of interest; FND must not be substituted.

## Statistical interpretation

Analyze the panels separately. For IC50, compare calculated energy with pIC50;
do not turn IC50 into a thermodynamic free energy. For the ITC panel, compare
against either pKd or the recorded standard-state `experimental_dg_kcal_mol`.
Because a more favorable calculated energy is more negative, the expected sign
is negative versus pIC50/pKd and positive versus experimental delta-G.

Spearman rank correlation is the primary metric. Pearson correlation and
leave-one-out results should be secondary because each panel has only seven
observations.

## Retrieval contract

- Target: the 14 exact BRD4-BD1 PDB-ligand pairs in `manifest.csv`.
- Organism: *Homo sapiens*, UniProt O60885, first bromodomain.
- Scope: bounded targeted retrieval, not an exhaustive PDB search.
- Primary database: RCSB Protein Data Bank.
- Endpoints: Data API `/rest/v1/core/entry/{PDB}`, coordinate downloads
  `/download/{PDB}.pdb` and `.cif`, and CCD `/ligands/download/{CCD}_ideal.sdf`.
- Local filters: exact manifest PDB and CCD; first coordinate model; highest
  heavy-atom-count/occupancy instance of the requested CCD; protein chains within
  8 A of that instance.
- Expected count: 14 structures, reported explicitly by the downloader.
- Activity access date: 2026-10-31; activity values are fixed in the manifest and
  do not come from potentially heterogeneous RCSB/BindingDB annotations.
