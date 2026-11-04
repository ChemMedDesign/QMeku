#!/usr/bin/env python3
"""Reproducible correlation analysis for the neutral BRD4 crystal benchmark."""

from __future__ import annotations

import csv
import itertools
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qmeku-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "results.csv"
OUT = ROOT / "analysis"
SEED = 20261031
BOOTSTRAPS = 20_000
PERMUTATIONS = 100_000
PANELS = ["pfi1_series", "wpf_shelf_series", "high_res_neutral", "all"]
METHODS = ["xtb", "mopac"]


def load_results() -> list[dict[str, object]]:
    with INPUT.open(newline="") as handle:
        rows: list[dict[str, object]] = list(csv.DictReader(handle))
    numeric = {
        "ic50_geomean_nM", "pIC50_exact", "resolution_A", "ligand_charge",
        "ligand_heavy_atoms", "pocket_residues", "pocket_atoms", "pocket_charge",
        "xtb_e_bind_kcal_mol", "mopac_e_bind_kcal_mol", "qmeku_seconds",
    }
    for row in rows:
        for field in numeric:
            row[field] = float(row[field])
    return rows


def subset(rows: list[dict[str, object]], panel: str) -> list[dict[str, object]]:
    return rows if panel == "all" else [row for row in rows if row["set"] == panel]


def correlation(kind: str, x: np.ndarray, y: np.ndarray) -> float:
    if kind == "spearman":
        return float(stats.spearmanr(x, y).statistic)
    return float(stats.pearsonr(x, y).statistic)


def bootstrap_ci(kind: str, x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float, int]:
    rng = np.random.default_rng(seed)
    values: list[float] = []
    n = len(x)
    for _ in range(BOOTSTRAPS):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 2 or np.unique(y[idx]).size < 2:
            continue
        value = correlation(kind, x[idx], y[idx])
        if math.isfinite(value):
            values.append(value)
    if not values:
        return math.nan, math.nan, 0
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high), len(values)


def spearman_permutation_p(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, str, int]:
    observed = abs(correlation("spearman", x, y))
    ranked_x = stats.rankdata(x)
    ranked_y = stats.rankdata(y)
    zx = (ranked_x - ranked_x.mean()) / ranked_x.std(ddof=0)
    zy = (ranked_y - ranked_y.mean()) / ranked_y.std(ddof=0)
    n = len(x)
    if n <= 8:
        extreme = total = 0
        for perm in itertools.permutations(range(n)):
            value = abs(float(np.mean(zx * zy[list(perm)])))
            extreme += value >= observed - 1e-12
            total += 1
        return extreme / total, "exact", total

    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(PERMUTATIONS):
        value = abs(float(np.mean(zx * zy[rng.permutation(n)])))
        extreme += value >= observed - 1e-12
    return (extreme + 1) / (PERMUTATIONS + 1), "monte_carlo", PERMUTATIONS


def holm_adjust(pvalues: list[float]) -> list[float]:
    order = np.argsort(pvalues)
    adjusted = np.empty(len(pvalues), dtype=float)
    running = 0.0
    m = len(pvalues)
    for rank, index in enumerate(order):
        running = max(running, (m - rank) * pvalues[index])
        adjusted[index] = min(1.0, running)
    return adjusted.tolist()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = load_results()
    OUT.mkdir(exist_ok=True)
    metrics: list[dict[str, object]] = []
    descriptives: list[dict[str, object]] = []
    size_rows: list[dict[str, object]] = []

    for panel_index, panel in enumerate(PANELS):
        data = subset(rows, panel)
        y = np.array([row["pIC50_exact"] for row in data], dtype=float)
        heavy = np.array([row["ligand_heavy_atoms"] for row in data], dtype=float)
        descriptives.append({
            "panel": panel,
            "n": len(data),
            "pIC50_min": y.min(),
            "pIC50_max": y.max(),
            "pIC50_span": y.max() - y.min(),
            "heavy_atoms_min": int(heavy.min()),
            "heavy_atoms_max": int(heavy.max()),
        })
        for method_index, method in enumerate(METHODS):
            x = np.array([row[f"{method}_e_bind_kcal_mol"] for row in data], dtype=float)
            descriptives[-1][f"{method}_energy_min"] = x.min()
            descriptives[-1][f"{method}_energy_max"] = x.max()
            for kind in ("spearman", "pearson"):
                estimate = correlation(kind, x, y)
                low, high, valid_boot = bootstrap_ci(
                    kind, x, y, SEED + panel_index * 100 + method_index * 10 + (kind == "pearson")
                )
                if kind == "spearman":
                    pvalue, p_method, n_perm = spearman_permutation_p(
                        x, y, SEED + panel_index * 10 + method_index
                    )
                else:
                    pvalue = float(stats.pearsonr(x, y).pvalue)
                    p_method, n_perm = "parametric_t", 0
                loo = []
                if len(x) >= 4:
                    for omitted in range(len(x)):
                        keep = np.arange(len(x)) != omitted
                        loo.append(correlation(kind, x[keep], y[keep]))
                residual_p = math.nan
                if kind == "pearson" and len(x) >= 4:
                    fit = np.polyfit(x, y, 1)
                    residual_p = float(stats.shapiro(y - np.polyval(fit, x)).pvalue)
                metrics.append({
                    "panel": panel,
                    "method": method,
                    "correlation": kind,
                    "n": len(x),
                    "estimate": estimate,
                    "ci95_boot_low": low,
                    "ci95_boot_high": high,
                    "valid_bootstraps": valid_boot,
                    "p_value_two_sided": pvalue,
                    "p_value_method": p_method,
                    "permutations": n_perm,
                    "p_holm_spearman_family": "",
                    "loo_min": min(loo) if loo else "",
                    "loo_max": max(loo) if loo else "",
                    "shapiro_residual_p": residual_p if math.isfinite(residual_p) else "",
                })

            energy_size = stats.spearmanr(x, heavy)
            efficiency_activity = stats.spearmanr(x / heavy, y)
            size_rows.append({
                "panel": panel,
                "method": method,
                "n": len(x),
                "rho_energy_vs_heavy_atoms": float(energy_size.statistic),
                "p_energy_vs_heavy_atoms": float(energy_size.pvalue),
                "rho_energy_per_heavy_atom_vs_pIC50": float(efficiency_activity.statistic),
                "p_energy_per_heavy_atom_vs_pIC50": float(efficiency_activity.pvalue),
            })

    spearman_indices = [i for i, row in enumerate(metrics) if row["correlation"] == "spearman"]
    adjusted = holm_adjust([float(metrics[i]["p_value_two_sided"]) for i in spearman_indices])
    for index, value in zip(spearman_indices, adjusted):
        metrics[index]["p_holm_spearman_family"] = value

    write_csv(OUT / "correlation_metrics.csv", metrics)
    write_csv(OUT / "descriptives.csv", descriptives)
    write_csv(OUT / "size_diagnostics.csv", size_rows)

    colors = {"pfi1_series": "#0072B2", "wpf_shelf_series": "#009E73", "high_res_neutral": "#D55E00"}
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True)
    for col, panel in enumerate(PANELS):
        data = subset(rows, panel)
        y = np.array([row["pIC50_exact"] for row in data], dtype=float)
        for row_index, method in enumerate(METHODS):
            ax = axes[row_index, col]
            x = np.array([row[f"{method}_e_bind_kcal_mol"] for row in data], dtype=float)
            for item, x_value, y_value in zip(data, x, y):
                ax.scatter(x_value, y_value, s=40, color=colors[str(item["set"])], edgecolor="black", linewidth=0.4)
                ax.annotate(str(item["pdb_id"]), (x_value, y_value), xytext=(3, 3), textcoords="offset points", fontsize=7)
            if len(x) >= 3 and np.unique(x).size > 1:
                line_x = np.linspace(x.min(), x.max(), 100)
                ax.plot(line_x, np.polyval(np.polyfit(x, y, 1), line_x), color="#555555", lw=1, ls="--")
            rho = correlation("spearman", x, y)
            r = correlation("pearson", x, y)
            ax.text(0.03, 0.04, f"n={len(x)}  ρ={rho:.2f}  r={r:.2f}", transform=ax.transAxes, fontsize=8)
            ax.set_title(panel.replace("_", " "))
            ax.set_xlabel(f"QMEku {method.upper()} E_bind (kcal/mol)")
            ax.set_ylabel("pIC50 experimental")
            ax.grid(alpha=0.2)
    fig.suptitle("BRD4-BD1 neutral crystal benchmark: calculated energy vs experimental activity", fontsize=14)
    fig.savefig(OUT / "correlation_panels.png", dpi=220)
    fig.savefig(OUT / "correlation_panels.pdf")
    plt.close(fig)

    print(f"Analyzed {len(rows)} complexes; outputs written to {OUT}")


if __name__ == "__main__":
    main()
