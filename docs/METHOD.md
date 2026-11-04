# The QMeku method

Everything is a single-point calculation with xtb or MOPAC, on the geometry of the pose, **without optimising it**. Energies in kcal/mol.

## 1. Pocket model (`qmeku/structure.py`)

- **Selection:** the protein residues with at least one atom (hydrogens included) closer than `--cutoff` Å to any atom of the ligand. With several poses and `--pocket common`, the union of the residues selected in any pose, so that all the poses use the same pocket and their energies are comparable.
- **Closure of the cut peptide bonds** with hydrogen link atoms, without optimising their position:
  - from the preceding residue (not included) the C and O are kept, and the Cα is replaced by an H along the C–Cα bond at 1.09 Å;
  - from the following one the N and H are kept, and the Cα (and the Cδ in proline) is replaced by an H along the N–C bond at 1.01 Å;
  - a cut disulfide bridge (precomputed symmetrically between `CYX`/`CYS` cysteines with no `HG` thiol closer than 2.5 Å) is closed with an H at 1.34 Å on the S.
- **Bridging residue:** a lone residue between two selected ones within the same chain is always included. If it were left out, its two neighbours would place a link atom on the same Cα and the two would end up 0.75 Å apart (observed in 3PTB with Ser217). In multi-chain receptors the labels carry a chain prefix (`A:ASP189`) to avoid collisions in oligomers.
- **Net charge and integrity checks:** the charge is counted by residue name (ASP/GLU/CYM −1; LYS/ARG/HIP +1) and by terminal ends verified by the absence of a peptide bond (protonated N terminus +1, including an N-terminal proline with two hydrogens on N; unprotonated OXT of the C terminus −1); the closures are neutral. The receptor must come with AMBER names (HID/HIE/HIP, CYX). The program rejects generic HIS, `CYX` without a disulfide partner, internal residues with spurious terminal atoms (`OXT`, `H1`–`H3`), atomic overlaps (< 0.5 Å) and any residue or subsystem with an odd number of electrons (even with `--total-only`).

## 2. Energies

For one method (xtb or MOPAC) and one phase (implicit solvent or gas), with `E(X)` the energy of a system `X`:

- **E_bind** (solvent) = E(P+L) − E(P) − E(L), with P the pocket and L the ligand.
- **E_gas** = the same difference in the gas phase. **ΔSolv** = E_bind − E_gas. That is, `E_bind = E_gas + ΔSolv`.
- **Per residue (`pair`)** = E(L+F) − E(L) − E(F), with F the residue closed on both sides. Also in the gas phase (`pair_gas`) and `dsolv = pair − pair_gas`.

xtb: GFN2-xTB with the ALPB solvent (water). MOPAC: PM6-D3H4X (dispersion and hydrogen-bond corrections) with COSMO (ε = 78.4); the energy is the heat of formation, and the complex − parts difference is the interaction energy because the atomisation references cancel.

### Concap correction (MFCC)

Two neighbouring fragments, F(i) and F(i+1), both contain the peptide group C(i)=O(i)–N(i+1)H(i+1) (each as the closure of the other). Adding their energies counts that group twice. Once per peptide bond between residues of the pocket, the interaction of the ligand with that group is subtracted (a formamide-like "concap", with hydrogen link atoms in place of the Cα atoms), split in halves between the two residues. This is the usual correction of molecular fragmentation with conjugate caps.

Effect measured in 3PTB (4 Å pocket, 10 peptide bonds between selected residues):

| | xtb uncorrected | xtb with concaps | MOPAC uncorrected | MOPAC with concaps |
|---|---|---|---|---|
| Σ residues (solvent) | −44.67 | −30.50 | +4.70 | +2.93 |
| E_bind | −25.36 | −25.36 | −17.41 | −17.41 |
| non-additive remainder (solvent) | +19.30 | +5.14 | −22.11 | −20.34 |
| remainder in the gas phase | +45.22 | +20.86 | +36.98 | +14.65 |

`--no-concap` switches it off. With xtb it corrects most of the non-additivity; with MOPAC a large remainder is left, originating in the COSMO solvation of each fragment.

### Balance: what the "non-additive remainder" is

`remainder = E_bind − Σ pair`. With the quantities of the log it decomposes into
`remainder = (E_gas − Σ pair_gas) + (ΔSolv − Σ dsolv)`:

- **Gas phase**: cooperative effects between residues (polarisation, charge transfer), which a calculation on pairs cannot capture.
- **Solvation**: each pair is computed with its own implicit solvent, so it carries a desolvation penalty of its own; the sum of those penalties exceeds that of the whole system (xtb, 3PTB: +97 over the pairs with concaps against +81 for the assembly).

The two parts have opposite signs in 3PTB (xtb: +20.9 in the gas phase and −15.7 from solvation). Saying that the remainder "is the desolvation" is therefore not correct: the pairs already include a desolvation, and an overestimated one.

## 3. Terms of the GFN2-xTB energy

**With xtb**, the following are read from the GFN2 energy summary: the electrostatic term (`es` = isotropic + anisotropic), the anisotropic XC (`xc`), the D4 dispersion (`disp`), the repulsion (`rep`), the ALPB solvation (`gsolv`) and a `rest` obtained by difference with the total (the remaining electronic energy), so that they **add up to the total**. For the interaction energy and for each residue they are subtracted in the same way as the total energy.

- They are terms of the GFN2 energy expression: analogous to, not identical with, classical vdW and Coulomb terms; they are not unique and they depend on the method.
- In the gas phase the interaction in 3PTB is almost all electrostatics and dispersion (xtb: ES −87.0, dispersion −19.6, repulsion +0.2, XC −0.2, rest −0.1; total −106.65).

**MOPAC** gives only the heat of formation and the COSMO dielectric energy; it does not separate electrostatics, dispersion or repulsion.

## 4. Convergence of xtb

Some systems (observed for an isolated ion pair, Asp194 + ligand, in the gas phase) do not converge in the SCF with the default parameters. `run_xtb` repeats the calculation with `--cycles 2000` and an electronic temperature of 600, 1000 and 2000 K, in that order, and records the change (`scf_fallbacks`). In a test with a system that converges either way, raising the temperature to 600 K changed the gas-phase energy by about 0.2 kcal/mol and left the ALPB one unchanged.

## 5. Dependence on the size of the pocket (3PTB, crystallographic pose, without concaps)

| Cutoff (Å) | Residues | Atoms | xtb/ALPB | PM6-D3H4X/COSMO |
|---|---|---|---|---|
| 3.0 | 11 | 149 | −26.5 | −18.3 |
| 4.0 | 16 | 225 | −25.4 | −17.4 |
| 5.0 | 24 | 348 | −24.9 | −22.3 |
| 6.0 | 26 | 389 | −25.5 | −22.5 |
| 7.0 | 33 | 515 | −26.2 | −20.3 |
| 8.0 | 50 | 752 | −25.1 | −20.5 |

xtb stays at −25.5 ± 1 from 3 Å onwards; MOPAC moves by about 5 kcal/mol between 4 and 5 Å and then lies between −20 and −22.5. The times of one complete interaction calculation (six single points) were 1 s with 149 atoms, 10 s with 348 and 72 s with 752.

Above roughly 500 atoms PM6-D3H4X becomes impracticable in this implementation: a 483-atom complex takes about 36 s, whereas a 609-atom complex (2478 electrons) did not complete a single SCF cycle in an hour of CPU. GFN2-xTB shows no such ceiling at these sizes.

## 6. What has not been validated

- The single-point energies have not been compared with a higher-level quantum reference method (DFT/ab initio) and, in the exploratory evaluation of 12 crystallographic complexes (`benchmarks/affinity_crystal/`, with neutral, anionic and cationic ligands), they do not predict experimental affinities across different proteins. Which of the two methods comes closer to a high-level reference is not known.
- The effect of relaxing the hydrogens has not been measured, nor that of a different scheme for closing the cuts.
- The per-residue decomposition has only been studied in detail on the development control (3PTB), and it measures the interaction of the ligand with each residue separately; it is not an exact partition of E_bind.
