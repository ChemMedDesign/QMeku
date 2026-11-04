# Analysis of the cases with a positive MOPAC energy

Date: 2026-11-01

## Main conclusion

The five positive cases share an energetic pattern, but **not** a single pattern of charge, size, target or geometric clash:

> With MOPAC/PM6-D3H4X, the stabilisation of the complex in the gas phase is insufficient to compensate the COSMO desolvation penalty. In four cases the gas-phase interaction is favourable but COSMO reverses the sign; in 4NUD the interaction is already unfavourable in the gas phase and the desolvation makes it worse.

There were no SCF failures, no odd-electron models and no heavy-atom contacts closer than 2.0 Å. The positive values are therefore reproducible results of the model and of the preparation used, not errors in reading the file.

## Energy balance

| PDB | q lig./pocket | MOPAC gas | ΔSolv COSMO | MOPAC final | xTB gas | ΔSolv ALPB | xTB final |
|---|---:|---:|---:|---:|---:|---:|---:|
| 3COW | 0 / −3 | −28.677 | **+75.407** | **+46.730** | −39.150 | +28.625 | −10.524 |
| 1AMW | −3 / 0 | −233.554 | **+257.325** | **+23.771** | −207.454 | +174.725 | −32.729 |
| 4NUD | 0 / 0 | **+4.451** | +18.735 | **+23.186** | −9.342 | +6.973 | −2.369 |
| 2V00 | 0 / −4 | −65.678 | **+82.865** | **+17.187** | −48.832 | +29.475 | −19.357 |
| 3MXF | +1 / −1 | −56.048 | **+68.670** | **+12.622** | −61.131 | +34.928 | −26.204 |

In the four cases with `E_gas < 0`, the ratio `ΔSolv/|E_gas|` with MOPAC was greater than one: 1AMW 1.10, 2V00 1.26, 3COW 2.63 and 3MXF 1.23. In five matched negative controls the ratio lay between 0.62 and 0.88; the gas-phase attraction kept dominating.

MOPAC imposed a larger desolvation penalty than xTB in all five positives. The median difference `ΔSolv(MOPAC) − ΔSolv(xTB)` was about 46.8 kcal/mol, against 21.9 kcal/mol in the controls. This identifies the COSMO/PM6-D3H4X treatment as the main common origin of the positive signs.

## Assessment of alternative patterns

### Patterns ruled out or weak

- **Ligand charge:** three ligands are neutral, one is −3 and one is +1. There is no common charge.
- **Pocket charge:** it ranges from 0 to −4. There are negative controls with equal or more extreme charges.
- **Size:** the positive pockets contain 254–483 atoms, within the range of the negative ones (212–577). The number of heavy atoms of the ligand does not separate the groups either.
- **Clashes:** the minimum heavy-atom distance was 2.62–3.15 Å and no positive case had heavy-atom contacts closer than 2.0 Å.
- **Target:** HSP90, endothiapepsin, pantothenate synthetase and BRD4 appear. Two are BRD4, but many BRD4 complexes came out negative.
- **Polarity:** TPSA was descriptively larger in the positives (median 130 against 87.7 Å²), driven by ADP/1AMW and 52H/3COW. The pattern is not universal: JQ1/3MXF is lipophilic (`cLogP ≈ 4.95`). It should be treated as a hypothesis, not a conclusion.
- **Excluded waters:** the positives had 7–16 crystallographic waters closer than 4 Å (median 10), but the controls showed 2–15. Removing waters may contribute case by case, although it does not by itself explain the selection.

Among the 46 original results, every positive showed a `MOPAC − xTB` discrepancy of at least +25.555 kcal/mol; the largest among the non-positives was +18.504 kcal/mol. This threshold is a diagnostic indicator of disagreement between methods within this dataset, not a universal rule.

## Case-by-case analysis

### 1AMW — HSP90/ADP

- ADP has charge −3, TPSA ≈ 241 Å², 14 acceptors and `cLogP ≈ −3.64`.
- The MOPAC gas-phase interaction is extraordinarily favourable (−233.554), but COSMO applies +257.325 kcal/mol.
- It is the clearest example of an extreme penalty associated with a polyanionic, heavily hydrated ligand.
- Interpretation: **strong evidence** that the positive sign comes from the COSMO desolvation, not from a bad pose.

### 2V00 — endothiapepsin/V15

- Neutral ligand in a pocket of charge −4.
- MOPAC goes from −65.678 in the gas phase to +17.187 through a penalty of +82.865.
- The crystal has seven waters closer than 4 Å; glycerol and acetates near the ligand were also excluded.
- Interpretation: **strong evidence** of a reversal by desolvation; **moderate evidence** of sensitivity to the explicit composition of the site.

### 3COW — pantothenate synthetase/52H

- It is the largest positive value: +46.730 kcal/mol.
- The neutral ligand is highly polar: TPSA ≈ 212 Å², 5 donors and 12 acceptors.
- COSMO contributes +75.407 against a gas-phase interaction of −28.677.
- The preparation records a geometric chain discontinuity `B:72 → B:81`; GLN71 is part of the immediate environment of the ligand, so the discontinuity is close to the site.
- Interpretation: **strong evidence** of excess desolvation and a **possibly important contribution** of the missing segment / chain break.

### 3MXF — BRD4/JQ1

- JQ1 has charge +1 and the pocket −1; the gas-phase interaction is favourable (−56.048).
- COSMO adds +68.670 and turns the result into +12.622.
- The polarity is not extreme (TPSA ≈ 71 Å²; `cLogP ≈ 4.95`), so the pattern cannot be reduced to "very polar ligands".
- There are 15 crystallographic waters closer than 4 Å. The crystallographic entry had three heavy atoms missing in LYS100; `pdb4amber --add-missing-atoms` rebuilt the whole side chain. It did not produce a direct clash, but its orientation is modelled and deserves a sensitivity test.
- Interpretation: **strong evidence** of a sensitivity specific to the PM6-D3H4X/COSMO balance; water and the side-chain repair are secondary hypotheses.

### 4NUD — BRD4/MS436

- It is the only case whose MOPAC interaction is already unfavourable in the gas phase: +4.451 kcal/mol.
- The desolvation adds +18.735, giving +23.186.
- It has no clashes. A neutral BRD4 control with the same charge and exactly 12 residues / 254 pocket atoms (4HBV) gives MOPAC gas −31.061, ΔSolv +19.268 and a total of −11.793.
- This rules out the sign being due to the size or charge of the model alone. The crystallographic entry had incomplete side chains, including LYS98 near the binding region; `pdb4amber --add-missing-atoms` rebuilt its atoms. The rebuilt orientation remains a local uncertainty.
- Interpretation: **strong evidence** of an unfavourable PM6 interaction for this microgeometry; protonation and tautomerism, waters and the repair of LYS98 call for a specific sensitivity test.

## Prioritised recommendations

1. Keep the five cases in the benchmark; do not truncate the values to zero nor remove them retrospectively.
2. Always report `E_gas` and `ΔSolv`. Flag `ΔSolv/|E_gas| ≥ 1` and `MOPAC − xTB > 20 kcal/mol` automatically for review.
3. Re-prepare 3COW, repairing or justifying the B72–B80 segment, and repeat with 4, 5 and 6 Å pockets.
4. Evaluate alternative microstates in ADP/1AMW, MS436/4NUD and JQ1/3MXF; keep the state fixed when comparing methods.
5. Test the sensitivity to the conserved crystallographic waters and other nearby heteroatoms, especially for 2V00, 3MXF and 4NUD, if the QMeku workflow comes to admit those components explicitly.
6. Repair and locally optimise nearby hydrogens and side chains before repeating 3MXF and 4NUD.
7. For publication, describe these points as **method-specific discrepancies dominated by the solvation balance**, not as evidence of an absence of experimental binding.

## Reproducible files

- `case_diagnostics.csv`: separated energies, properties, contacts, waters and provenance.
- `<PDB>/diagnostics/gas_split/`: the complete diagnostic calculation of each positive case.
- `analysis/controls/<PDB>/gas_split/`: five matched negative controls.
- `analyze_cases.py`: recompiles the structural diagnostics.

The critical-assessment framework followed Kassis T., Agarwal V., He Y., Patel D. and Brueckner A. M. (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, arXiv:2609.00065, DOI: [10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065).
