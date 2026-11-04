# QMeku

**Which residues actually hold your ligand — and how much of that the solvent takes back.**

![version](https://img.shields.io/badge/version-1.0-blue)
![python](https://img.shields.io/badge/python-%E2%89%A5%203.9-blue)
![licence](https://img.shields.io/badge/licence-MIT-green)
![tests](https://img.shields.io/badge/tests-32-informational)

QMeku takes a receptor and a set of docking poses and computes the ligand–pocket interaction energy with
**semiempirical quantum mechanics** — GFN2-xTB through [xtb](https://github.com/grimme-lab/xtb) and
PM6-D3H4X through [MOPAC](https://openmopac.net/) — instead of a fitted scoring function. It then splits that
energy **residue by residue**, separates the **gas-phase attraction from the desolvation penalty**, and writes
a PyMOL scene with the pocket coloured by each residue's contribution.

One command, two engines, no parameters to fit:

```bash
qmeku examples/3ptb/receptor.pdb examples/3ptb/crystal.sdf --outdir out
```

```
E_bind  GFN2-xTB (ALPB): -25.36   PM6-D3H4X (COSMO): -17.41
Balance:  E_int gas -106.65 / -115.05   dSolv +81.29 / +97.64
Residue   d(A)  heavy atoms      GFN2-xTB  PM6-D3H4X
ASP189    1.95  N2-OD2 (2.87)      -13.14     -6.28
```

Trypsin and benzamidine, crystallographic pose: the salt bridge to **ASP189** carries half the interaction,
and the solvent gives back +81 kcal/mol of it. That is the kind of statement QMeku is for.

### Why you might want it

- **A decomposition that adds up.** The per-residue table closes against E_bind through two explicit rows,
  `Σ residues` and `non-additive remainder` — nothing is hidden in a fitted constant. The MFCC/concap
  correction keeps the shared backbone from being counted twice (in 3PTB it cuts the remainder from +19.3 to
  +5.1 kcal/mol).
- **Two independent methods on the same pocket.** Agreement between GFN2-xTB and PM6-D3H4X is information;
  so is disagreement, and QMeku makes both visible instead of averaging them away.
- **The gas/solvation balance, in the open.** `E_bind = E_gas + ΔSolv`. For charged ligands these are two
  huge numbers that nearly cancel; QMeku prints both, so you can see when a result rests on that cancellation.
- **Fast enough to use.** The example above: **3 s** and **114 MB** of peak memory on 8 cores, full
  decomposition included.
- **Output you can look at.** CSV files for analysis, a readable log, and `pymol/` scenes with the receptor
  coloured by residue contribution.
- **46 benchmark calculations on 44 crystal structures ship with the code**, results and failures included.

### What it is not

> **The energies are not validated against experimental affinity.** They are single-point calculations on a
> fixed geometry: no optimisation, no entropy, no ligand or protein reorganisation, no explicit waters. Across
> heterogeneous targets they do **not** correlate with experimental free energy. Use them to compare poses and
> to read which residues matter — not to predict a Kd.

That boundary is documented rather than glossed over: see [Benchmarks](#benchmarks) and
[Limitations](#limitations), and `benchmarks/BENCHMARK_REPORT.md` for the full report, including the five
complexes where PM6-D3H4X comes out repulsive and why.

---

## Minimum requirements

**Required**

| | Minimum | Notes |
|:--|:--|:--|
| Python | 3.9 | The package itself uses only the standard library plus RDKit |
| [RDKit](https://www.rdkit.org/) | any recent release | Ligand perception, hydrogens, formal charges, SDF input. Tested with 2026.03.6 |
| A quantum engine | **xtb or MOPAC — at least one** | Both is better: you get both methods on the same pocket. Tested with xtb 6.6.1 and 6.7.1, MOPAC 23.2.2 and 23.2.5 |
| OS | 64-bit Linux | Developed and tested on Linux x86-64. Nothing in the code is platform-specific and both engines exist for macOS, but macOS and Windows are untested |
| Hardware | 1 core, ~0.5 GB free RAM, 20 MB disk | Measured on the bundled example (one pose, xtb, full per-residue decomposition): 3 s and 114 MB peak on 8 cores. More cores mean more poses and residues in parallel (`--workers`) |

**Optional**

| | For what |
|:--|:--|
| [PDB2PQR](https://github.com/Electrostatics/pdb2pqr) | `--prepare`: protonate the receptor at pH 7.4 with AMBER names |
| [PyMOL](https://pymol.org/) | Rendering `.png`/`.pse` images from the scene scripts. The `.pml` scripts are written without it |
| numpy, scipy, pandas, matplotlib, Pillow | Re-running the statistical analysis of `benchmarks/` (`pip install "qmeku[benchmarks]"`) |
| AmberTools (`pdb4amber`, `reduce`) | Re-preparing the benchmark receptors from the crystals |

**Pocket size is the real limit, not the hardware.** PM6-D3H4X becomes impracticable somewhere above 500
atoms: a 483-atom complex takes about 36 s, while a 609-atom complex did not finish a single SCF cycle in an
hour of CPU. GFN2-xTB is unaffected at these sizes. Keep MOPAC pockets under ~500 atoms, or lower `--cutoff`.

## Installation

### 1. The engines (install at least one)

```bash
# conda / mamba, the easiest route for both
conda install -c conda-forge xtb        # GFN2-xTB, LGPL-3.0
conda install -c conda-forge mopac      # PM6-D3H4X, academic licence — read it at openmopac.net
conda install -c conda-forge rdkit      # or: pip install rdkit
```

QMeku finds them on the `PATH`, or through the variables `XTB_BIN`, `MOPAC_BIN` and `PYMOL_BIN` if you keep
them somewhere else. Nothing is bundled or redistributed: see `NOTICE.md`.

### 2. QMeku

```bash
git clone https://github.com/ChemMedDesign/qmeku.git
cd qmeku
pip install .                      # creates the `qmeku` command
```

Variants:

```bash
pip install -e .                   # editable, for development
pip install ".[benchmarks]"        # also numpy/scipy/pandas/matplotlib/Pillow for the benchmark analysis
bin/qmeku --help                   # no installation at all: run it from the source tree
```

### 3. Check that it works

```bash
qmeku --version                                   # qmeku 1.0
qmeku --examples                                  # usage examples, no calculation
make test PYTHON=/path/to/python-with-rdkit       # 32 unit tests
qmeku examples/3ptb/receptor.pdb examples/3ptb/crystal.sdf --outdir out --method xtb
```

The last command is the real check: a few seconds, and `out/` should contain `qmeku.log`, `poses.csv`,
`residues.csv` and `pymol/`. Tests that need xtb or MOPAC skip themselves when the program is missing, so a
green test run does **not** prove the engines are installed — the example does.

## Usage

```
qmeku RECEPTOR.pdb POSES.sdf --outdir RESULTS
```

- **RECEPTOR.pdb** — with hydrogens and AMBER residue names (HID/HIE/HIP instead of HIS, CYX for disulfides).
  Amino acids only: no water, ions or cofactors. `--prepare` protonates it for you with PDB2PQR.
- **POSES.sdf** — one molecule per pose, in 3D, placed in the site, with correct bonds and formal charges.
  Each conformation of a molecule counts as one pose. The ligand charge comes from RDKit unless you set
  `--ligand-charge`.

```bash
# one pose, xtb and MOPAC, with the gas/solvation balance
qmeku examples/3ptb/receptor.pdb examples/3ptb/crystal.sdf --outdir out

# tens of poses: xtb only, no gas phase (faster), per-residue average over the poses with E_bind <= 0
qmeku receptor.pdb poses.sdf --outdir out --method xtb --no-gas --summary-max-ebind 0

# log in Spanish
qmeku receptor.pdb poses.sdf --outdir out --lang es
```

| Option | Meaning |
|---|---|
| `--method xtb\|mopac\|both` | Semiempirical program (default `both`; if one is missing, the program says so and continues with the other) |
| `--solvent water\|none` | Implicit solvent (water by default) or gas phase only |
| `--cutoff X` | Pocket cutoff in Å (4.0 by default) |
| `--pocket common\|per-pose` | `common` (the default with several poses): the same residues in all of them, which is needed to compare energies between poses. `per-pose`: each pose gets its own |
| `--no-gas` | Skips the gas phase (faster; no gas/solvation balance) |
| `--examples` | Prints usage examples and exits |
| `--no-concap` | Does not subtract the peptide groups shared by neighbouring fragments (see `docs/METHOD.md`) |
| `--summary-max-ebind X` | Excludes from the per-residue average the poses whose E_bind (first method) is larger than X |
| `--ligand-charge N` | Net charge of the ligand (by default, the one from RDKit) |
| `--prepare` | Protonates the receptor with PDB2PQR (AMBER names, pH 7.4) before running |
| `--total-only` | E_bind only, without the per-residue decomposition |
| `--no-pymol` | Does not write the PyMOL view (`OUTDIR/pymol/`), which is generated by default together with the per-residue decomposition (see `docs/OUTPUTS.md`) |
| `--pymol-poses N` | Number of best poses drawn in PyMOL (3 by default; 0 = all) |
| `--pymol-method xtb\|mopac\|both` | Method(s) drawn in PyMOL, among those computed (default `both`: one scene per method) |
| `--top N` | Residues shown in the tables of the log (0 = all) |
| `--workers N` | Simultaneous calculations (all cores by default) |
| `--timeout N` | Maximum time in seconds per xtb/MOPAC/pdb2pqr calculation (600 by default; 0 = no limit) |
| `--lang en\|es` | Language of the log (`en` by default) |
| `--quiet` | Does not print the log to the screen (it is always written to `qmeku.log`) |

## How it works

For each pose:

1. **Pocket model.** The residues with at least one atom closer than the cutoff (4.0 Å by default) to the
   ligand. Cut peptide bonds are closed with hydrogen link atoms; bridging residues and odd-electron models
   are handled explicitly (`docs/METHOD.md`).
2. **Energies**, with xtb (GFN2-xTB, ALPB water) and/or MOPAC (PM6-D3H4X, COSMO):
   - **E_bind** = E(pocket + ligand) − E(pocket) − E(ligand), in implicit solvent;
   - the same in the **gas phase** plus the **solvation change**: `E_bind = E_gas + ΔSolv`;
   - the **contribution of each residue** — the ligand's interaction with that residue alone, again with its
     own gas phase and solvation change;
   - with xtb, the **GFN2 energy terms** (electrostatics, D4 dispersion, repulsion, anisotropic XC,
     solvation) for the total and for each residue.
3. **Outputs**: CSV files, a readable log, and the PyMOL scenes.

## Outputs

In `RESULTS/`:

| File | Content |
|---|---|
| `qmeku.log` | Readable report: E_bind per pose and method, gas/solvation balance, GFN2-xTB terms, per-residue table whose rows add up to E_bind, summary over all the poses |
| `poses.csv` | One row per pose: E_bind, E_gas, ΔSolv, sum over residues, non-additive remainder, xtb terms, convergence retries |
| `residues.csv` | One row per pose and residue: contribution, gas phase, ΔSolv, xtb terms, distances and closest atoms |
| `residues_summary.csv` | With several poses: per residue, contacts, mean, standard deviation and median |
| `pymol/` | One `.pml` scene per pose and method, with the receptor coloured by residue contribution |

Every column is described in `docs/OUTPUTS.md`.

## How to read the results

**E_bind** is the ligand–pocket interaction in implicit solvent. It serves to **compare** poses that share the
same pocket; it is not an affinity. For two poses to be comparable they must use the same residues
(`--pocket common`, the default with several poses).

**The balance matters.** `E_bind = E_gas + ΔSolv`. In the gas phase the interaction of a cationic ligand with
a carboxylate is enormous (≈ −107 kcal/mol in 3PTB with xtb) and the solvent compensates it with a large
desolvation penalty (≈ +81). E_bind is the difference between two large numbers, and that is worth knowing
before trusting a small one.

**The sum over residues is not E_bind.** The table includes the rows `Σ residues`, `non-additive remainder`
and `E_bind`, which add up. The remainder gathers two effects of opposite sign, not only the desolvation:
cooperative effects between residues (gas phase) and solvation (each pair is computed with its own implicit
solvent and carries a desolvation penalty which, taken separately, is overestimated). The concap correction
also avoids counting twice the backbone shared by neighbouring residues.

**GFN2-xTB terms.** Their sum is E_bind. They are terms of the GFN2 energy expression — not unique
thermodynamic partitions, and not equivalents of classical force-field terms. MOPAC does not separate them.

## Benchmarks

`benchmarks/` ships four sets — **46 computed complexes over 44 distinct crystal structures**, two of which
appear in two sets for a sensitivity comparison — with the prepared inputs, every result table, the analysis
scripts and a README per set. `benchmarks/BENCHMARK_REPORT.md` is the full report.

| Set | What it tests | Headline result |
|:--|:--|:--|
| `affinity_crystal` | 12 heterogeneous complexes, Kd/Ki | No correlation with experimental free energy across different targets |
| `brd4_neutral_series` | 17 BRD4-BD1 complexes, neutral ligands | PM6-D3H4X ranks them: Spearman ρ = −0.723, 95% CI [−0.876, −0.382], permutation *p* = 0.0014 (Holm-adjusted 0.011). GFN2-xTB is inconclusive: ρ = −0.397, *p* = 0.12 |
| `brd4_expanded` | ITC/IC50 series, chain and protonation sensitivity | How much the answer moves when you change a defensible preparation choice |
| `mopac_positive_cases` | The 5 complexes where PM6-D3H4X comes out **repulsive** | Traced to the COSMO desolvation penalty exceeding the gas-phase attraction — not to execution errors |

Two things to take from this. The ranking result holds **within one target and one assay**, and does not
transfer to a heterogeneous set. And the five positive cases are shipped, not dropped: a method's failures
are part of its description.

The crystals themselves are not redistributed. `download_crystals.py` in each set re-fetches them from the
RCSB, and the SHA-256 of every file is recorded in the corresponding `metadata.json`; the `prepared/` inputs
that the calculations actually used **are** included, so every number here is reproducible byte for byte
without any download.

## Limitations

- **No predictive validation of affinity.** Single-point energies do not correlate with experimental free
  energy across heterogeneous targets. Which of the two methods is closer to a higher-level quantum reference
  is not known: they agree on the dominant residue and differ in magnitude. MOPAC depends on pocket size more
  than xtb does (in 3PTB, −17.4 with 16 residues and −22.95 with 24; xtb −25.4 and −25.2).
- **Unoptimised geometry.** Receptor and ligand hydrogens are not relaxed; the effect has not been measured.
- **Simple closure of the cuts** (hydrogen link atoms), and a pocket made of amino acids only — no waters,
  ions or cofactors.
- **No entropy and no reorganisation energy.**
- **Clashing poses** give enormous, meaningless energies; there is no clash flag in the semiempirical part
  (use `--summary-max-ebind`).
- **No per-atom decomposition of the ligand**; for each residue only the closest atoms are reported.
- **The xtb SCF sometimes fails to converge** with the default parameters (seen for an isolated ion pair in
  the gas phase). It is retried at a higher electronic temperature, and this is reported in the log and in
  `poses.csv` (`scf_fallbacks`).
- **Pocket size ceiling with MOPAC**, as described under Minimum requirements.

## Layout

```
bin/qmeku            launcher (no installation needed)
qmeku/               Python package: cli, ligand (RDKit and PDB), structure (pocket, caps, concaps),
                     engines (xtb, MOPAC), analysis (pose), report (CSV and log), pymol_view
examples/3ptb/       receptor, crystallographic ligand, SDF of one pose and of 20 synthetic poses
benchmarks/          four benchmark sets, each with its own README, and BENCHMARK_REPORT.md
docs/METHOD.md       method, formulas, decisions and validation
docs/OUTPUTS.md      columns of each file
tests/               32 unit tests
NOTICE.md            provenance and licences
```

## Citing QMeku

The article describing QMeku and these benchmarks is in preparation. Until it is published, please cite the
repository:

> SB2CS lab. *QMeku 1.0: approximate binding energy and per-residue decomposition of docking poses with
> semiempirical methods.* 2026. https://github.com/ChemMedDesign/qmeku

Please also cite the engines you used: GFN2-xTB and ALPB for xtb, PM6-D3H4X and COSMO for MOPAC. The
references are in `docs/METHOD.md`.

## Licence

MIT, © 2026 SB2CS lab — see `LICENSE`. QMeku invokes xtb, MOPAC, RDKit and PDB2PQR but redistributes none of
them; their licences and the provenance of the code are in `NOTICE.md`. **MOPAC's licence is academic: check
that your use is covered before you rely on PM6-D3H4X results.**

Developed by the **SB2CS lab**. Bug reports and questions are welcome through the issue tracker.
