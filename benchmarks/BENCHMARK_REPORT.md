# Complete report of the crystallographic benchmarking of QMeku

**Status of the document:** local, reproducible technical report; it is not a manuscript approved for submission.
**Closing date:** 2026-10-31 (local time in Chile; 2026-11-01 UTC).
**Software evaluated:** QMeku 1.0, the local code current at the time of the calculations.
**Scope:** from the first heterogeneous affinity benchmark to the BRD4-BD1 extension and the sensitivity analyses for chain selection and protonation.

---

## 1. Aim

The aim was to assess whether the interaction energy computed by QMeku correlates with experimental measurements of affinity or potency for crystallographic protein–ligand complexes. The aim was **not** to reproduce the experimental free energy numerically, but to determine whether the method preserves the relative order of affinity within comparable sets.

QMeku computes a single-point interaction energy:

\[
E_{\mathrm{bind}} = E(\mathrm{pocket+ligand}) - E(\mathrm{pocket}) - E(\mathrm{ligand})
\]

This quantity is not a standard binding free energy. It does not explicitly include entropy, ligand reorganisation, full protein reorganisation, conformational populations or the contribution of discrete structural waters.

## 2. Common computational protocol

The benchmarks were run with the following protocol:

- crystallographic pose of the ligand, with no optimisation of the heavy atoms;
- protein preparation with AmberTools `pdb4amber --reduce`;
- ligand preparation with RDKit and an explicit charge/protonation decision;
- pocket formed by the residues with at least one atom closer than 4.0 Å to the ligand;
- GFN2-xTB with ALPB water;
- PM6-D3H4X with COSMO water through MOPAC;
- `--total-only --no-gas`;
- no per-residue decomposition in the main calculation;
- two internal processes per complex;
- preservation of the crystallographic heavy-atom coordinates, with a heavy-atom RMSD of 0.000 Å in the validated preparations.

The prepared receptors contain amino acids, hydrogens and AMBER names, with no waters, ions or cofactors. The ligand charges come from the documented chemical state or from an explicit decision for pH 7.4. The workflow and its limitations are described in [affinity_crystal/README.md](affinity_crystal/README.md) and [brd4_expanded/README.md](brd4_expanded/README.md).

---

## 3. First benchmark: heterogeneous crystallographic set

### 3.1. Composition

The first benchmark gathered twelve cases: six central `Kd` complexes, three stress cases — two `Ki` and one `Kd` —, three BRD4 complexes used to assess ranking, and one development control, 3PTB. Some complexes belong simultaneously to an experimental category and to an interpretative subset; the correlations were therefore computed by strata and not through a single confirmatory pool.

The structures were downloaded from the RCSB PDB, keeping the PDB, mmCIF, metadata, the CCD template of the ligand and the crystallographic coordinates. The selection and the experimental labels are recorded in [affinity_crystal/manifest.csv](affinity_crystal/manifest.csv). [L1]

### 3.2. Recalculation after changes in the code

After QMeku was modified, the twelve complexes were recomputed with local version 1.0 and the protocol above.

During the review it was found that 1BCU contained PDB insertion codes that produced spurious internal termini during the preparation. After the numbering and the geometric breaks were corrected:

- the pocket changed from 253 atoms/charge −6 to 248 atoms/charge −5;
- xTB changed from −61.550 to −21.941 kcal/mol;
- MOPAC changed from −8.376 to −11.950 kcal/mol.

The incident is documented in [affinity_crystal/README.md](affinity_crystal/README.md). [L5]

### 3.3. Current results of the first benchmark

**Table 1.** Consolidated values after the recalculation. The numerical source is [affinity_crystal/results.csv](affinity_crystal/results.csv). [L2]

| PDB | CCD ligand | Experimental label | Experimental ΔG (kcal/mol) | xTB (kcal/mol) | MOPAC (kcal/mol) | Main use |
|---|---|---:|---:|---:|---:|---|
| 2V00 | V15 | Kd 220 µM | −4.99 | −19.357 | +17.187 | central Kd |
| 1ADL | ACD | Kd 4.4 µM | −7.31 | −6.486 | −0.724 | central Kd |
| 3P5O | EAM | Kd 55.2 nM | −9.90 | −26.259 | −19.688 | central Kd / BRD4 |
| 3COW | 52H | Kd 125 nM | −9.42 | −10.524 | +46.730 | central Kd; two mutations |
| 1AMW | ADP | Kd 29 µM | −6.19 | −32.729 | +23.771 | central Kd; ligand −3 |
| 1BM7 | FLF | Kd 30 nM | −10.26 | −11.517 | −6.893 | central Kd; oligomeric interface |
| 5F63 | 5W2 | Kd 12 nM | −10.80 | −40.934 | −30.168 | BRD4 |
| 4WIV | 3P2 | Kd 550 nM | −8.54 | −23.839 | −23.357 | BRD4 |
| 1E66 | HUX | Ki 0.13 nM | −13.48 | −40.028 | −21.524 | stress / Ki |
| 2HB1 | 512 | Ki 160 µM | −5.18 | −39.504 | −36.819 | stress / Ki; ligand −2 |
| 1BCU | PRL | Kd 530 µM | −4.47 | −21.941 | −11.950 | stress; insertion-code correction |
| 3PTB | BEN | Ki 18.4 µM | −6.46 | −25.362 | −17.412 | development control |

### 3.4. Correlation of the first benchmark

The correlations were analysed against the `experimental ΔG`; in this convention a **positive** correlation is expected, because a more favourable experimental binding and a more favourable computed interaction are both more negative.

**Table 2.** Pre-stratified correlations. Source: [affinity_crystal/analysis/correlation_metrics.csv](affinity_crystal/analysis/correlation_metrics.csv). [L3]

| Panel | n | Method | Pearson | Spearman | Interpretation |
|---|---:|---|---:|---:|---|
| central Kd | 6 | xTB | −0.293 | −0.200 | no signal; very wide bootstrap 95% CI |
| central Kd | 6 | MOPAC | +0.314 | +0.543 | uncertain, not significant |
| all Kd | 9 | xTB | +0.103 | +0.200 | no useful correlation |
| all Kd | 9 | MOPAC | +0.287 | +0.417 | weak and uncertain |
| initial BRD4 | 3 | xTB | +0.873 | +1.000 | correct order, descriptive only |
| initial BRD4 | 3 | MOPAC | +0.546 | +0.500 | partial order, descriptive only |
| full pool | 12 | xTB | +0.175 | +0.259 | descriptive; does not validate absolute affinity |
| full pool | 12 | MOPAC | +0.155 | +0.259 | descriptive; does not validate absolute affinity |

For the `Kd` panels, 20,000 bootstrap resamples and seed 20261031 were used. The intervals were wide and the conclusions sensitive to the exclusion of observations. The manifest of the analysis is in [affinity_crystal/analysis/analysis_manifest.json](affinity_crystal/analysis/analysis_manifest.json). [L4]

### 3.5. Initial BRD4 subset

The three initial BRD4 structures have the following experimental order:

| PDB | Kd | pKd | xTB | MOPAC |
|---|---:|---:|---:|---:|
| 5F63 | 12 nM | 7.921 | −40.934 | −30.168 |
| 3P5O | 55.2 nM | 7.258 | −26.259 | −19.688 |
| 4WIV | 550 nM | 6.260 | −23.839 | −23.357 |

xTB reproduces the experimental order of these three cases exactly. The size `n=3` prevents this observation from being treated as inferential validation. These complexes are considered an external panel because their `Kd` values do not come from the same ITC study used later.

### 3.6. Exploratory comparison with IC50

The extraction of exact PDB–ligand annotations from BindingDB gave coverage for only 3 of the 12 initial cases: 1ADL, 1BCU and 5F63. The comparison was made through pIC50, without treating IC50 as a thermodynamic constant.

| Method | n | Pearson with pIC50 | Spearman with pIC50 |
|---|---:|---:|---:|
| xTB | 3 | −0.827 | −0.500 |
| MOPAC | 3 | −0.868 | −0.500 |

Although the sign of Pearson is the expected one — a more negative energy against a higher pIC50 — `n=3` and the heterogeneity of the assays make the result purely descriptive. Sources: [affinity_crystal/analysis_ic50/ic50_comparison.csv](affinity_crystal/analysis_ic50/ic50_comparison.csv) and [affinity_crystal/analysis_ic50/ic50_metrics.csv](affinity_crystal/analysis_ic50/ic50_metrics.csv). [L6]

---

## 4. Extended BRD4-BD1 benchmark

### 4.1. Design

An independent folder, `brd4_expanded`, was built with two panels of seven complexes each:

1. **IC50 panel:** the same chemical series, human BRD4-BD1 and the Oxford B AlphaScreen assay with four replicates, taken from Wang et al., DOI [10.1021/acschembio.7b00638](https://doi.org/10.1021/acschembio.7b00638).
2. **Kd/ITC panel:** direct affinities for BRD4-BD1 measured by ITC in Karim et al., DOI [10.1021/acs.jmedchem.1c01096](https://doi.org/10.1021/acs.jmedchem.1c01096).

The activities, errors, CCD codes, charge decisions and provenance are kept in [brd4_expanded/manifest.csv](brd4_expanded/manifest.csv). [L7]

### 4.2. Chemical preparation

The seven ligands of the IC50 panel contain a basic tertiary amine and were represented with charge +1 at pH 7.4. The seven ligands of the ITC panel also ended up with charge +1.

JQ1/3MXF had initially been marked as neutral. The preparation failed visibly because the official CCD `JQ1` carries a formal charge of +1, consistent with the diazepin-10-ium name. The manifest was corrected to keep the CCD protomer directly, without adding another protonation site.

In 6CJ1, `F2S` was selected specifically, as the ligand of interest. The `FND` stereoisomer present in the same entry was not used.

The pipeline completed 14/14 cases.

### 4.3. Results of the IC50 panel

For IC50 the following was used:

\[
pIC_{50}=-\log_{10}(IC_{50}[\mathrm{M}])
\]

IC50 was not converted into `ΔG`, because IC50 depends on the design of the assay and is not necessarily an equilibrium constant.

**Table 3.** Primary IC50 panel. Source: [brd4_expanded/results.csv](brd4_expanded/results.csv). [L9]

| PDB | Ligand | IC50 (µM) | pIC50 | xTB | MOPAC |
|---|---|---:|---:|---:|---:|
| 6CIS | JWG-047 | 0.30 ± 0.05 | 6.523 | −22.957 | −14.930 |
| 6CIY | JWG-069 | 0.41 ± 0.01 | 6.387 | −20.061 | −16.519 |
| 6CD5 | XMD17-26 | 0.75 ± 0.11 | 6.125 | −18.990 | −16.638 |
| 5WA5 | XMD11-50/LRRK2-IN-1 | 1.04 ± 0.18 | 5.983 | −14.001 | −9.174 |
| 5W55 | JWG-048 | 2.41 ± 0.23 | 5.618 | −23.370 | −20.491 |
| 6CD4 | JWG-046 | 5.14 ± 0.59 | 5.289 | −16.295 | −14.085 |
| 6CJ1 | JWG-071/F2S | 6.31 ± 1.93 | 5.200 | −21.910 | −14.300 |

### 4.4. Results of the Kd/ITC panel

For Kd the following were used:

\[
pK_d=-\log_{10}(K_d[\mathrm{M}])
\]

and

\[
\Delta G^{\circ}_{\mathrm{exp}}=RT\ln K_d
\]

with `T = 298.15 K` for the tabulated values.

**Table 4.** Primary Kd/ITC panel before the multi-chain corrections. Source: [brd4_expanded/results.csv](brd4_expanded/results.csv). [L9]

| PDB | Ligand | Kd (nM) | pKd | exp. ΔG | xTB | MOPAC |
|---|---|---:|---:|---:|---:|---:|
| 3MXF | (+)-JQ1 | 70.9 ± 6.6 | 7.149 | −9.753 | −26.204 | +12.622 |
| 4OGI | BI2536 | 112 ± 5 | 6.951 | −9.483 | −30.567 | −14.993 |
| 5VBO | LRRK2-IN-1 | 183 ± 14 | 6.738 | −9.192 | −21.061 | −17.111 |
| 7K6G | ERK5-IN-1 | 213.2 ± 10.8 | 6.671 | −9.101 | −19.747 | −16.796 |
| 5VBP | Ro3280 | 266.7 ± 23.8 | 6.574 | −8.969 | −33.604 | −15.523 |
| 5V67 | Volasertib | 317 ± 17 | 6.499 | −8.866 | −32.142 | −17.976 |
| 5LRQ | XMD8-92 | 735 ± 38 | 6.134 | −8.368 | −24.199 | −7.654 |

### 4.5. Primary correlations of the extended panel

The following correlations use pActivity; a **negative** sign is therefore expected: a more negative energy should accompany a higher pActivity.

**Table 5.** Primary correlations. Source: [brd4_expanded/correlation_metrics.csv](brd4_expanded/correlation_metrics.csv). [L10]

| Panel | Method | n | Pearson | Pearson *p* | Spearman | Spearman *p* |
|---|---|---:|---:|---:|---:|---:|
| IC50 | xTB | 7 | −0.097 | 0.836 | −0.107 | 0.819 |
| IC50 | MOPAC | 7 | −0.010 | 0.982 | −0.286 | 0.535 |
| Kd/ITC | xTB | 7 | −0.019 | 0.968 | +0.143 | 0.760 |
| Kd/ITC | MOPAC | 7 | +0.447 | 0.315 | +0.321 | 0.482 |

No statistically significant correlation and no consistent ranking were found on the complete panels.

### 4.6. Sensitivity through 3MXF/JQ1

3MXF produced a sign discrepancy between the methods:

- xTB: −26.204 kcal/mol;
- MOPAC: +12.622 kcal/mol.

The +1 charge is not a reading error: it comes from the `JQ1` CCD entry and was kept deliberately. The difference between the methods nevertheless identifies 3MXF as a case sensitive to the electronic representation, the implicit solvent and the parameterisation.

Excluding 3MXF as a sensitivity analysis only:

| Method | n | Pearson with pKd | Spearman with pKd |
|---|---:|---:|---:|
| xTB | 6 | −0.065 | +0.257 |
| MOPAC | 6 | −0.684 | −0.086 |

The exclusion does not produce a robust improvement of the ranking. 3MXF must not be excluded from the primary analysis merely because its result is unfavourable.

---

## 5. Investigation of the extreme values 5VBP and 5V67

### 5.1. Decomposition of the xTB terms

The xTB values of 5VBP and 5V67 were among the most negative of the panel, even though their Kd values were not the best.

| xTB term (kcal/mol) | 5VBP | 5V67 |
|---|---:|---:|
| Electrostatics | −27.355 | −9.253 |
| D4 dispersion | −35.285 | −31.236 |
| Anisotropic XC | −0.786 | −0.967 |
| Repulsion | +0.152 | +0.133 |
| ALPB solvation | +15.233 | −4.389 |
| Remaining electronic | +14.437 | +13.570 |
| **Total E_bind** | **−33.604** | **−32.142** |

In 5VBP, electrostatics and dispersion dominate, partly offset by desolvation. In 5V67, dispersion and a favourable solvation term dominate. These terms belong to the GFN2-xTB energy expression and must not be read as unique thermodynamic partitions.

### 5.2. Multi-chain artefact in 4OGI and 5VBP

The initial downloader selected every protein chain lying closer than 8 Å to the ligand. In 4OGI and 5VBP this included chains A and B, although the selected ligand was associated with chain B. The original preparation therefore incorporated additional crystal-packing contacts.

A sensitivity analysis was built that keeps only the protein chain matching the chain identifier of the ligand. The original benchmark was left unoverwritten.

**Table 6.** Effect of the chain selection. Source: [brd4_expanded/chain_matched_sensitivity/comparison.csv](brd4_expanded/chain_matched_sensitivity/comparison.csv). [L11]

| PDB | Model | Pocket residues/atoms | Pocket charge | xTB | MOPAC |
|---|---|---:|---:|---:|---:|
| 4OGI | chains A+B | 20 / 400 | −1 | −30.567 | −14.993 |
| 4OGI | chain B only | 15 / 296 | 0 | **−21.927** | **−20.061** |
| 5VBP | chains A+B | 18 / 380 | 0 | −33.604 | −15.523 |
| 5VBP | chain B only | 13 / 276 | +1 | **−23.836** | **−19.371** |

The correction makes xTB 8.640 kcal/mol less favourable in 4OGI and 9.768 kcal/mol less favourable in 5VBP. This confirms that the extreme xTB value of 5VBP was strongly influenced by the inclusion of the second chain.

MOPAC changes in the opposite direction: −5.068 kcal/mol in 4OGI and −3.848 kcal/mol in 5VBP. This does not contradict the structural correction, because `E_bind` is a non-additive difference between three calculations with implicit solvent; removing atoms does not guarantee a monotonic change.

### 5.3. Protonation microstate of 5V67

Volasertib contains a disubstituted piperazine with two protonatable nitrogens. The primary calculation protonated the internal atom N39. A second protomer with the same total charge of +1 was generated, protonating N40. The heavy atoms kept an RMSD of 0.000 Å.

**Table 7.** Sensitivity to the protonation microstate. Source: [brd4_expanded/protonation_sensitivity/5V67_alternate_piperazine_N/comparison.json](brd4_expanded/protonation_sensitivity/5V67_alternate_piperazine_N/comparison.json). [L12]

| 5V67 protomer | xTB | MOPAC |
|---|---:|---:|
| N39 protonated, primary | −32.142 | −17.976 |
| N40 protonated, alternative | **−28.880** | **−18.687** |
| Alternative − primary | +3.262 | −0.711 |

The microstate explains part of the xTB magnitude, but not the whole favourable value. The size of Volasertib — 45 heavy atoms, about 620 Da — and the dispersion remain important factors. These properties were obtained from the prepared SDF structure used in the calculation. [L13]

### 5.4. Correlation after the sensitivity analyses

Replacing 4OGI and 5VBP by the single-chain calculations, keeping the primary protomer of 5V67:

| Corrected panel | Method | n | Pearson | Spearman |
|---|---|---:|---:|---:|
| full Kd/ITC | xTB | 7 | +0.151 | +0.250 |
| full Kd/ITC | MOPAC | 7 | +0.374 | +0.036 |

The complete panel still shows no correlation in the expected direction.

Excluding 3MXF in addition, as a post hoc sensitivity:

| Corrected panel without 3MXF | Method | n | Pearson | *p* | Spearman | *p* |
|---|---|---:|---:|---:|---:|---:|
| Kd/ITC | xTB | 6 | +0.389 | 0.445 | +0.714 | 0.111 |
| Kd/ITC | MOPAC | 6 | **−0.850** | **0.032** | −0.543 | 0.266 |

With the alternative protomer of 5V67, MOPAC keeps a similar Pearson signal (`r = −0.829`, `p = 0.041`), while Spearman stays at −0.543. This divergence indicates that the signal depends on the magnitudes and does not yet constitute a robust ranking. It also rests on `n=6` and on a post hoc exclusion; it must not be presented as a confirmatory result.

---

## 6. Integrated interpretation

1. **QMeku does not reproduce absolute affinities.** This is consistent with its computing a single-point electronic interaction, not a free energy.
2. **The heterogeneity between targets dominates the first benchmark.** Mixing proteins, charges, sizes and assay types removes any simple relation between `E_bind` and the experimental `ΔG`.
3. **The first BRD4 subset showed an encouraging ranking, but with only three points.** xTB ordered 5F63, 3P5O and 4WIV correctly, without sufficient statistical power.
4. **The homogeneous extension did not confirm that signal.** The panels of seven IC50 and seven Kd/ITC complexes show no significant correlation in their primary analysis.
5. **The definition of the receptor is material.** Including a crystal-packing chain changed xTB by about 9–10 kcal/mol in 4OGI and 5VBP.
6. **The protonation is material.** Changing the protonated nitrogen of Volasertib modified xTB by 3.262 kcal/mol.
7. **xTB and MOPAC respond differently.** 3MXF even changes sign between methods, and the chain corrections move xTB and MOPAC in different directions.
8. **The pockets are not identical across crystals.** The 4 Å criterion yields different numbers of residues, atoms and charges; this limits direct comparison even within BRD4.

---

## 7. Recommendations for the next publishable benchmark

### 7.1. Structures and receptor

- use the single-chain receptors for 4OGI and 5VBP;
- audit automatically whether the ligand is in contact with a second crystallographic chain;
- align the BRD4-BD1 domains and define a common set of pocket residues for the comparative analysis;
- keep an analysis with an individual pocket as a sensitivity test as well;
- do not mix packing contacts with the biological site unless there is explicit structural justification.

### 7.2. Ligands

- compute or document the dominant protonation microstates at the experimental pH;
- for 5V67, keep both protomers until a better-supported pKa/microstate assignment is available;
- keep 3MXF in the primary analysis and treat its exclusion only as a sensitivity test;
- keep the stereochemistry, the heavy-atom coordinates and the CCD–PDB trail.

### 7.3. Statistics

- use Spearman as the primary ranking metric;
- use Pearson as the secondary metric of linear relation;
- report bootstrap confidence intervals and leave-one-out analyses;
- analyse IC50 and Kd separately;
- do not convert IC50 into `ΔG`;
- keep the three initial BRD4 cases as external validation, not as part of the fit of the ITC panel;
- do not choose exclusions on the basis of whether they improve the correlation.

### 7.4. Development of the method

- assess the energy per heavy atom, or the interaction efficiency, only as a secondary analysis, not as an automatic substitute for `E_bind`;
- study the effect of a common pocket, of structural waters and of histidine microstates;
- compare with a higher-level reference method on a small subset;
- document the effects of size, charge, dispersion and implicit solvent separately;
- keep the original results alongside every correction, to avoid retrospective bias.

---

## 8. Limitations

- The sample sizes are small: `n=3`, `n=6` or `n=7` in the main analyses.
- The structures are single points and the added hydrogens were not relaxed by a full optimisation.
- The solvent is implicit and does not represent discrete structural waters.
- The energies depend on the size and composition of the selected pocket.
- The sum of terms, or of per-residue contributions, is not necessarily equal to `E_bind`, because of non-additive effects.
- The initial Kd values and the extended Kd/ITC values do not all come from the same study.
- The corrected correlations without 3MXF are exploratory and post hoc.
- A null correlation does not show that the method is useless for other tasks, such as discriminating poses under one and the same receptor and pocket.

---

## 9. File layout and reproducibility

### Initial benchmark

- experimental selection: [affinity_crystal/manifest.csv](affinity_crystal/manifest.csv)
- consolidated results: [affinity_crystal/results.csv](affinity_crystal/results.csv)
- correlations: [affinity_crystal/analysis/correlation_metrics.csv](affinity_crystal/analysis/correlation_metrics.csv)
- IC50 analysis: [affinity_crystal/analysis_ic50/ic50_metrics.csv](affinity_crystal/analysis_ic50/ic50_metrics.csv)
- protocol and decisions: [affinity_crystal/README.md](affinity_crystal/README.md)

### Extended BRD4 benchmark

- selection and activities: [brd4_expanded/manifest.csv](brd4_expanded/manifest.csv)
- consolidated results: [brd4_expanded/results.csv](brd4_expanded/results.csv)
- primary correlations: [brd4_expanded/correlation_metrics.csv](brd4_expanded/correlation_metrics.csv)
- multi-chain sensitivity: [brd4_expanded/chain_matched_sensitivity/comparison.csv](brd4_expanded/chain_matched_sensitivity/comparison.csv)
- protonation sensitivity: [brd4_expanded/protonation_sensitivity/5V67_alternate_piperazine_N/comparison.json](brd4_expanded/protonation_sensitivity/5V67_alternate_piperazine_N/comparison.json)
- protocol: [brd4_expanded/README.md](brd4_expanded/README.md)

Main commands:

```bash
cd /path/to/qmeku

# Initial benchmark
python benchmarks/affinity_crystal/run_benchmark.py
python benchmarks/affinity_crystal/collect_results.py
python benchmarks/affinity_crystal/analyze_correlations.py

# Extended BRD4
python benchmarks/brd4_expanded/run_all.py --method both --workers 2
python benchmarks/brd4_expanded/analyze_results.py

# Later sensitivity analyses
python benchmarks/brd4_expanded/run_chain_matched_sensitivity.py
python benchmarks/brd4_expanded/run_5v67_protonation_sensitivity.py
```

The `raw/` directories of the downloaded crystals are not distributed with this repository. They are recovered with the `download_crystals.py` script of each set, and the SHA-256 of every downloaded file is recorded in the corresponding `metadata.json`.

---

## 10. Local provenance record

| ID | Local evidence | Content used |
|---|---|---|
| L1 | `affinity_crystal/manifest.csv` | PDB codes, ligands, affinity, paths and initial classification |
| L2 | `affinity_crystal/results.csv` | Energies and pocket properties of the initial benchmark |
| L3 | `affinity_crystal/analysis/correlation_metrics.csv` | Initial correlations, bootstrap and leave-one-out |
| L4 | `affinity_crystal/analysis/analysis_manifest.json` | Number of rows, seed, resamples and analytical versions |
| L5 | `affinity_crystal/README.md` | Protocol, recalculation and correction of 1BCU |
| L6 | `affinity_crystal/analysis_ic50/*` | Coverage and exploratory IC50 correlations |
| L7 | `brd4_expanded/manifest.csv` | IC50/Kd activities, errors, charges, DOIs and pActivity |
| L9 | `brd4_expanded/results.csv` | Primary energies of the 14 BRD4 complexes |
| L10 | `brd4_expanded/correlation_metrics.csv` | Primary correlations and sensitivity without 3MXF |
| L11 | `brd4_expanded/chain_matched_sensitivity/comparison.csv` | Chain correction of 4OGI and 5VBP |
| L12 | `brd4_expanded/protonation_sensitivity/5V67_alternate_piperazine_N/comparison.json` | Protonation sensitivity of Volasertib |
| L13 | `brd4_expanded/kd_itc_series/5V67/prepared/ligand.sdf` | Prepared structure used to count heavy atoms and compute the molecular mass |

---

## 11. References

1. Wang et al. *Structural and Atropisomeric Factors Governing the Selectivity of Pyrimido-benzodiazipinones as Inhibitors of Kinases and Bromodomains*. ACS Chemical Biology. DOI: [10.1021/acschembio.7b00638](https://doi.org/10.1021/acschembio.7b00638).
2. Karim RM, Bikowitz MJ, Chan A, et al. *Differential BET Bromodomain Inhibition by Dihydropteridinone and Pyrimidodiazepinone Kinase Inhibitors*. Journal of Medicinal Chemistry. 2021;64:15772–15786. DOI: [10.1021/acs.jmedchem.1c01096](https://doi.org/10.1021/acs.jmedchem.1c01096).
3. RCSB Protein Data Bank. Structures, metadata and Chemical Component Dictionary used through the identifiers recorded in the manifests. [https://www.rcsb.org/](https://www.rcsb.org/).
4. Kassis T, Agarwal V, He Y, Patel D, Brueckner AM. *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065, current version consulted on 2026-10-31. DOI: [10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065).

---

## 12. Conclusion

The first benchmark showed that QMeku has no general correlation with absolute affinities across heterogeneous targets. The initial BRD4 subset gave a favourable ranking, especially with xTB, but it was limited to three observations. The extension to two homogeneous panels of seven complexes did not confirm a robust primary correlation.

The later revisions showed that apparently minor structural decisions can change `E_bind` by several kcal/mol: the inclusion of a second crystallographic chain altered xTB by about 9–10 kcal/mol in 4OGI and 5VBP, and the protonation microstate of Volasertib changed xTB by 3.262 kcal/mol. The priority methodological next step is therefore not to fit QMeku to the experimental values, but to standardise the receptor, the pocket and the protonation before assessing the ranking again.
