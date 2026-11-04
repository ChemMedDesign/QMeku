#!/usr/bin/env python3
"""Compare QMeku energies with exact PDB/ligand IC50 annotations."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qmeku-ic50")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats


BENCHMARK = Path(__file__).resolve().parent
RESULTS = BENCHMARK / "results.csv"
OUTPUT = BENCHMARK / "analysis_ic50"
TEMPERATURE_K = 298.15
R_KCAL_MOL_K = 1.98720425864083e-3
METHODS = {
    "xTB/ALPB": "xtb_e_bind_kcal_mol",
    "MOPAC/COSMO": "mopac_e_bind_kcal_mol",
}


def entry_json_by_pdb() -> dict[str, Path]:
    return {
        path.parents[1].name.upper(): path
        for path in BENCHMARK.glob("*/*/raw/entry.json")
    }


def extract_ic50_measurements(results: pd.DataFrame) -> pd.DataFrame:
    entries = entry_json_by_pdb()
    rows: list[dict[str, object]] = []
    for result in results.itertuples(index=False):
        path = entries.get(result.pdb_id.upper())
        if path is None:
            continue
        entry = json.loads(path.read_text(encoding="utf-8"))
        for value in entry.get("rcsb_binding_affinity", []):
            if value.get("type") != "IC50" or value.get("comp_id") != result.ligand_ccd:
                continue
            if value.get("unit") != "nM":
                raise ValueError(f"Unexpected IC50 unit in {path}: {value.get('unit')}")
            ic50_nm = float(value["value"])
            rows.append({
                "pdb_id": result.pdb_id,
                "ligand_ccd": result.ligand_ccd,
                "ic50_nM": ic50_nm,
                "pIC50": 9.0 - math.log10(ic50_nm),
                "reference_sequence_identity_percent": value.get("reference_sequence_identity"),
                "provenance": value.get("provenance_code"),
                "source_url": value.get("link"),
                "rcsb_entry_json": str(path.relative_to(BENCHMARK)),
            })
    return pd.DataFrame(rows)


def aggregate_ic50(measurements: pd.DataFrame, results: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (pdb_id, ligand_ccd), group in measurements.groupby(["pdb_id", "ligand_ccd"], sort=True):
        values = group["ic50_nM"].to_numpy(float)
        geometric_mean = float(np.exp(np.mean(np.log(values))))
        pic50 = 9.0 - math.log10(geometric_mean)
        # This is only an IC50-derived concentration scale. IC50 is not an
        # equilibrium constant, so the result must not be called experimental ΔG.
        dg_equivalent = R_KCAL_MOL_K * TEMPERATURE_K * math.log(geometric_mean * 1e-9)
        source = results.loc[results.pdb_id == pdb_id].iloc[0]
        rows.append({
            "pdb_id": pdb_id,
            "set": source["set"],
            "ligand_ccd": ligand_ccd,
            "ic50_n": len(values),
            "ic50_geomean_nM": geometric_mean,
            "ic50_min_nM": float(values.min()),
            "ic50_max_nM": float(values.max()),
            "pIC50": pic50,
            "ic50_dg_equivalent_kcal_mol": dg_equivalent,
            "reference_sequence_identity_min_percent": float(group["reference_sequence_identity_percent"].min()),
            "xtb_e_bind_kcal_mol": source["xtb_e_bind_kcal_mol"],
            "mopac_e_bind_kcal_mol": source["mopac_e_bind_kcal_mol"],
            "source_url": group["source_url"].iloc[0],
        })
    return pd.DataFrame(rows).sort_values("pdb_id").reset_index(drop=True)


def descriptive_metrics(comparison: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    # Higher pIC50 means stronger potency; more negative QMeku energy is
    # interpreted as stronger interaction. A useful association is negative.
    x = comparison["pIC50"].to_numpy(float)
    for method, column in METHODS.items():
        y = comparison[column].to_numpy(float)
        pearson = stats.pearsonr(x, y)
        spearman = stats.spearmanr(x, y)
        rows.extend([
            {
                "method": method,
                "n": len(comparison),
                "x_variable": "pIC50",
                "correlation": "Pearson",
                "estimate": pearson.statistic,
                "p_value": pearson.pvalue,
                "interpretation": "descriptive_only_n_lt_5",
            },
            {
                "method": method,
                "n": len(comparison),
                "x_variable": "pIC50",
                "correlation": "Spearman",
                "estimate": spearman.statistic,
                "p_value": spearman.pvalue,
                "interpretation": "descriptive_only_n_lt_5",
            },
        ])
    return pd.DataFrame(rows)


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    png = OUTPUT / f"{stem}.png"
    fig.savefig(png, dpi=300, facecolor="white", transparent=False)
    with Image.open(png) as image:
        image.convert("RGB").save(png, dpi=(300, 300))
    fig.savefig(OUTPUT / f"{stem}.pdf", facecolor="white", transparent=False)


def plot(comparison: pd.DataFrame, metrics: pd.DataFrame) -> None:
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.7), sharex=True, layout="constrained")
    colors = {"1ADL": "#0072B2", "5F63": "#D55E00", "1BCU": "#009E73"}
    for ax, (method, column), panel in zip(axes, METHODS.items(), ("A", "B")):
        for row in comparison.itertuples(index=False):
            ax.scatter(row.pIC50, getattr(row, column), s=60, color=colors[row.pdb_id],
                       edgecolor="black", linewidth=0.6, zorder=3)
            ax.annotate(row.pdb_id, (row.pIC50, getattr(row, column)),
                        xytext=(5, 4), textcoords="offset points", fontsize=8)
        r = metrics[(metrics.method == method) & (metrics.correlation == "Pearson")].iloc[0]
        rho = metrics[(metrics.method == method) & (metrics.correlation == "Spearman")].iloc[0]
        ax.text(0.04, 0.96, f"r = {r.estimate:.2f}\nρ = {rho.estimate:.2f}\nn = 3",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#777777"})
        ax.set_title(method, fontsize=10)
        ax.set_xlabel("pIC50 (geometric mean; molar)")
        ax.grid(True, color="#DDDDDD", linewidth=0.6, zorder=0)
        ax.text(-0.12, 1.04, panel, transform=ax.transAxes, fontsize=12, fontweight="bold")
    axes[0].set_ylabel("QMeku interaction energy (kcal/mol)")
    fig.suptitle("Exact PDB/ligand IC50 subset — exploratory only", fontsize=11)
    save_figure(fig, "ic50_comparison")
    plt.close(fig)


def write_report(results: pd.DataFrame, comparison: pd.DataFrame, metrics: pd.DataFrame) -> None:
    covered = set(comparison.pdb_id)
    missing = results.loc[~results.pdb_id.isin(covered), ["pdb_id", "ligand_ccd"]]
    lines = [
        "# QMeku comparison with IC50",
        "",
        "## Outcome",
        "",
        f"Exact PDB/ligand IC50 annotations were available for {len(comparison)} of {len(results)} complexes. "
        "This is too small for validation or inferential correlation; the coefficients below are descriptive diagnostics only.",
        "",
        "Replicate IC50 values were combined with a geometric mean because IC50 is analyzed on a logarithmic concentration scale. "
        "All individual values remain in `ic50_measurements.csv`.",
        "",
        "| PDB | Ligand | IC50 values (nM) | Geometric mean (nM) | pIC50 | xTB | MOPAC | Sequence identity |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    measurements = pd.read_csv(OUTPUT / "ic50_measurements.csv")
    for row in comparison.itertuples(index=False):
        vals = measurements.loc[measurements.pdb_id == row.pdb_id, "ic50_nM"]
        values = ", ".join(f"{value:g}" for value in vals)
        lines.append(
            f"| {row.pdb_id} | {row.ligand_ccd} | {values} | {row.ic50_geomean_nM:.3g} | "
            f"{row.pIC50:.3f} | {row.xtb_e_bind_kcal_mol:.3f} | {row.mopac_e_bind_kcal_mol:.3f} | "
            f"{row.reference_sequence_identity_min_percent:.0f}% |"
        )
    lines += ["", "## Descriptive correlations", ""]
    for method in METHODS:
        r = metrics[(metrics.method == method) & (metrics.correlation == "Pearson")].iloc[0]
        rho = metrics[(metrics.method == method) & (metrics.correlation == "Spearman")].iloc[0]
        lines.append(f"- {method}: Pearson r={r.estimate:.3f}; Spearman ρ={rho.estimate:.3f}; n=3.")
    lines += [
        "",
        "A favorable association is negative because stronger potency raises pIC50 while stronger computed interaction makes the energy more negative. "
        "With only three heterogeneous targets, neither the coefficient nor its p-value is suitable as publication evidence.",
        "",
        "## Complexes without exact IC50 annotation",
        "",
        ", ".join(f"{row.pdb_id}/{row.ligand_ccd}" for row in missing.itertuples(index=False)) + ".",
        "",
        "## Interpretation constraints",
        "",
        "- IC50 is assay-dependent and is not an equilibrium binding constant. pIC50 is the defensible comparison scale.",
        "- `ic50_dg_equivalent_kcal_mol` is supplied only as a plotting convenience (`RT ln(IC50 in M)`); it is not an experimental binding free energy.",
        "- Converting IC50 to Ki requires the assay mechanism and quantities such as substrate/tracer concentration and Km or Kd.",
        "- 1ADL has only 91% reference-sequence identity in the BindingDB annotation and should be sensitivity-flagged.",
        "- The three observations span unrelated targets and special structural contexts; no global calibration should be fitted.",
        "",
        "## Recommendation for a publication",
        "",
        "Present this IC50 analysis as a secondary, explicitly exploratory panel. For a defensible IC50 benchmark, expand the crystal set prospectively using one assay family per target, exact ligand identity, documented construct/sequence, and preferably at least 20–30 independent complexes. "
        "Keep Kd, Ki, and IC50 in separate strata rather than converting and pooling them silently.",
        "",
        "## Provenance",
        "",
        "IC50 values were read from each downloaded RCSB `entry.json` field `rcsb_binding_affinity`, whose provenance is BindingDB. "
        "Raw rows, including source URLs and sequence identity, are in `ic50_measurements.csv`.",
    ]
    (OUTPUT / "ic50_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(RESULTS)
    measurements = extract_ic50_measurements(results)
    if measurements.empty:
        raise RuntimeError("No exact PDB/ligand IC50 annotations found")
    measurements.to_csv(OUTPUT / "ic50_measurements.csv", index=False, float_format="%.8g")
    comparison = aggregate_ic50(measurements, results)
    comparison.to_csv(OUTPUT / "ic50_comparison.csv", index=False, float_format="%.8g")
    metrics = descriptive_metrics(comparison)
    metrics.to_csv(OUTPUT / "ic50_metrics.csv", index=False, float_format="%.8g")
    plot(comparison, metrics)
    write_report(results, comparison, metrics)
    print(f"IC50 coverage: {len(comparison)}/{len(results)} complexes")
    print(comparison[["pdb_id", "ic50_geomean_nM", "pIC50", *METHODS.values()]].to_string(index=False))


if __name__ == "__main__":
    main()
