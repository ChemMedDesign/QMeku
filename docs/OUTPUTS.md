# Output files

Energies in kcal/mol, distances in Å. An empty cell means that the datum does not exist (for example the `es`, `disp`… terms
with MOPAC, the gas-phase columns with `--no-gas`, or `d_solv` / `solv_*` when `--solvent none` is used). `{m}` is `xtb` or `mopac`.

## `poses.csv` (one row per pose)

| Column | Meaning |
|---|---|
| `pose` | Name of the molecule in the SDF (with `_c<n>` if it has several conformations) |
| `pocket_residues`, `pocket_atoms`, `pocket_charge` | Size and net charge of the pocket (bridging residues included) |
| `ligand_charge` | Net charge of the ligand (from RDKit or `--ligand-charge`) |
| `solvent`, `seconds` | Implicit solvent and calculation time |
| `scf_fallbacks` | xtb calculations repeated with a higher electronic temperature, e.g. `xtb gas FL ASP194 etemp 600` |
| `{m}_e_bind` | **Approximate binding energy** in implicit solvent |
| `{m}_e_gas`, `{m}_d_solv` | Gas-phase interaction and solvation change; `e_bind = e_gas + d_solv` |
| `{m}_pair_sum` | Sum of the per-residue contributions (with the concap correction) |
| `{m}_residual` | Non-additive remainder: `e_bind − pair_sum` |
| `{m}_pair_sum_raw`, `{m}_concap_sum` | Uncorrected sum and sum of the concap corrections; `pair_sum = pair_sum_raw − concap_sum` |
| `{m}_pair_sum_gas`, `{m}_residual_gas` | The same in the gas phase: `e_gas − pair_sum_gas` |
| `{m}_pair_dsolv_sum`, `{m}_residual_solv` | Sum of the per-residue solvation changes and its remainder: `d_solv − pair_dsolv_sum`; `residual = residual_gas + residual_solv` |
| `{m}_solv_{es,xc,disp,rep,gsolv,rest}` | GFN2 terms of `e_bind` (xtb only; `gsolv` is the ALPB solvation) |
| `{m}_gas_{es,xc,disp,rep,gsolv,rest}` | Terms of `e_gas` (xtb only; `gsolv` = 0) |

The xtb terms add up to the total: `es + xc + disp + rep + gsolv + rest = e_bind` (and `e_gas` in the gas phase).

## `residues.csv` (one row per pose and pocket residue)

| Column | Meaning |
|---|---|
| `residue`, `resname`, `resnum` | Label (`ASP189` with a single chain or `A:ASP189` in multi-chain receptors), name and number of the residue |
| `index` | Position of the residue in the file |
| `in_contact` | 1 if the closest distance is ≤ `--cutoff` |
| `min_dist`, `ligand_atom`, `ligand_atom_index`, `receptor_atom` | Closest distance between any pair of atoms (hydrogens included) and the atoms that define it |
| `min_dist_heavy`, `heavy_ligand_atom`, `heavy_receptor_atom` | The same between heavy atoms (the one quoted for hydrogen bonds and salt bridges) |
| `n_atoms_within_cutoff` | Atoms of the residue closer than the cutoff to the ligand |
| `{m}_pair` | **Contribution of the residue**: interaction of the ligand with the residue alone, in solvent, with the concap correction |
| `{m}_pair_raw` | The same without the concap correction |
| `{m}_pair_gas`, `{m}_dsolv` | In the gas phase and the solvation change: `pair = pair_gas + dsolv` |
| `{m}_solv_*`, `{m}_gas_*` | GFN2 terms of the residue (xtb only) |

With `--pocket common`, every pose has one row for each residue of the common pocket, including those far from that pose.

## `residues_summary.csv` (with several poses)

Per residue: `n_poses`, `n_contact` and `contact_fraction` (poses with the residue in contact), `min_dist_mean` and `min_dist_best`,
and per method `{m}_mean`, `{m}_sd` (population standard deviation, divisor `n_poses`) and `{m}_median` of the contribution over the poses. With `--summary-max-ebind X` only the poses with E_bind (first method) ≤ X are counted. A residue absent from a pose
(per-pose pocket) counts as 0 in the energies and does not enter `min_dist_mean`.

## `qmeku.log`

A readable report with a header (versions, methods, pocket, warning) and, per pose: the binding energy with each method, the
gas/solvation balance, the xtb terms in the gas phase, the GFN2 terms in solvent, the per-residue table
(contribution, distances and atoms) with the rows `Σ residues`, `non-additive remainder` and `E_bind`, and the
gas/solvation breakdown. At the end, the summary over all the poses and the mean contribution per residue.

## PyMOL representation (`pymol/`)

When the per-residue decomposition is computed (that is, without `--total-only`), qmeku writes by default into `OUTDIR/pymol/` the
ligand–receptor view coloured by `{m}_pair` (red = most negative, light grey ≈ 0, light blue = positive). It is switched off with
`--no-pymol`. The `--pymol-poses` best poses are drawn (lowest E_bind of the first method; 3 by default, 0 = all) for each
computed method (`--pymol-method xtb|mopac|both`, `both` by default), and with several poses one further scene with the mean contribution (`residues_summary.csv`).

| File | Content |
|---|---|
| `<pose>_<method>.pml` | Self-contained PyMOL script (plain commands, relative paths): `pymol <pose>_<method>.pml` from that folder |
| `<pose>_<method>.png` | Image; only if PyMOL is installed (`pymol` on the PATH or the `PYMOL_BIN` variable) |
| `summary_<method>.pml/.png` | The same with the mean over poses (several poses) |
| `receptor.pdb`, `ligand_<pose>.sdf` | Coordinates loaded by the `.pml` files |
| `<pose>_<method>_receptor_Bfactor.pdb` | Receptor with the energy of each residue in the B factor (only with PyMOL) |

The colour scale is common to the poses drawn for one method. The ligand is pale yellow (small spheres and
sticks), only the polar hydrogens of the binding site and of the ligand are shown, the most important residues are labelled
(`Asp189`) and the polar ligand–residue contacts (≤ 3.5 Å) are drawn. If PyMOL is not installed, only the `.pml` files are written.
