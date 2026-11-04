#!/usr/bin/env python3
"""Exploratory correlation analysis for the QMeku crystal benchmark."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qmeku-correlation")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats


BENCHMARK = Path(__file__).resolve().parent
INPUT = BENCHMARK / "results.csv"
OUTPUT = BENCHMARK / "analysis"
SEED = 20261031
N_BOOT = 20_000
METHODS = {
    "xTB/ALPB": "xtb_e_bind_kcal_mol",
    "MOPAC/COSMO": "mopac_e_bind_kcal_mol",
}


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    png = OUTPUT / f"{stem}.png"
    fig.savefig(png, dpi=300, facecolor="white", transparent=False)
    # Matplotlib PNGs commonly retain an opaque alpha channel. Deliver an
    # explicit RGB file so its background cannot change in another context.
    with Image.open(png) as image:
        image.convert("RGB").save(png, dpi=(300, 300))
    fig.savefig(OUTPUT / f"{stem}.pdf", facecolor="white", transparent=False)


def bootstrap_ci(x: np.ndarray, y: np.ndarray, kind: str) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = []
    for _ in range(N_BOOT):
        indices = rng.integers(0, len(x), len(x))
        bx, by = x[indices], y[indices]
        if len(np.unique(bx)) < 2 or len(np.unique(by)) < 2:
            continue
        value = stats.pearsonr(bx, by).statistic if kind == "pearson" else stats.spearmanr(bx, by).statistic
        if np.isfinite(value):
            estimates.append(value)
    if not estimates:
        return float("nan"), float("nan")
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def leave_one_out_range(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    values = [stats.pearsonr(np.delete(x, index), np.delete(y, index)).statistic for index in range(len(x))]
    return min(values), max(values)


def metric_rows(data: pd.DataFrame, group_name: str, inferential: bool) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    x = data["experimental_dg_kcal_mol"].to_numpy(float)
    for method, column in METHODS.items():
        y = data[column].to_numpy(float)
        pearson = stats.pearsonr(x, y)
        spearman = stats.spearmanr(x, y)
        slope, intercept = np.polyfit(x, y, 1)
        fitted = intercept + slope * x
        residuals = y - fitted
        shapiro = stats.shapiro(residuals) if len(data) >= 4 else None
        loo_min, loo_max = leave_one_out_range(x, y) if len(data) >= 4 else (float("nan"), float("nan"))
        for kind, result in (("Pearson", pearson), ("Spearman", spearman)):
            ci = bootstrap_ci(x, y, kind.lower()) if inferential else (float("nan"), float("nan"))
            rows.append({
                "group": group_name,
                "n": len(data),
                "method": method,
                "correlation": kind,
                "estimate": result.statistic,
                "p_value_two_sided": result.pvalue if inferential else float("nan"),
                "bootstrap_ci95_low": ci[0],
                "bootstrap_ci95_high": ci[1],
                "bootstrap_resamples": N_BOOT if inferential else 0,
                "linear_slope": slope,
                "linear_intercept": intercept,
                "residual_shapiro_w": shapiro.statistic if shapiro else float("nan"),
                "residual_shapiro_p": shapiro.pvalue if shapiro else float("nan"),
                "pearson_leave_one_out_min": loo_min,
                "pearson_leave_one_out_max": loo_max,
                "inferential_interpretation": inferential,
            })
    return rows


def scatter_figure(data: pd.DataFrame, metrics: pd.DataFrame) -> None:
    kd = data[data["affinity_type"] == "Kd"].copy()
    styles = {
        "core_absolute": ("#0072B2", "o", "Core Kd"),
        "brd4_rank": ("#D55E00", "s", "BRD4 rank"),
        "stress": ("#009E73", "^", "Stress Kd"),
    }
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.7), sharex=True, sharey=True, layout="constrained")
    for ax, (method, column), panel in zip(axes, METHODS.items(), ("A", "B")):
        for set_name, subset in kd.groupby("set", sort=False):
            color, marker, label = styles[set_name]
            ax.scatter(
                subset["experimental_dg_kcal_mol"], subset[column],
                s=52, c=color, marker=marker, edgecolors="black", linewidths=0.55,
                label=label, zorder=3,
            )
        x = kd["experimental_dg_kcal_mol"].to_numpy(float)
        y = kd[column].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(x.min(), x.max(), 100)
        ax.plot(grid, intercept + slope * grid, color="#333333", linestyle="--", linewidth=1.1, label="OLS, Kd")
        for _, row in kd.iterrows():
            ax.annotate(row["pdb_id"], (row["experimental_dg_kcal_mol"], row[column]),
                        xytext=(4, 4), textcoords="offset points", fontsize=7)
        m = metrics[(metrics["group"] == "Kd_all") & (metrics["method"] == method)]
        r = m[m["correlation"] == "Pearson"].iloc[0]
        rho = m[m["correlation"] == "Spearman"].iloc[0]
        ax.text(0.03, 0.97, f"r = {r['estimate']:.2f}\nρ = {rho['estimate']:.2f}\nn = {len(kd)}",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#777777", "alpha": 0.9})
        ax.set_title(method, fontsize=10)
        ax.set_xlabel("Experimental ΔG-equivalent (kcal/mol)")
        ax.grid(True, color="#DDDDDD", linewidth=0.6, zorder=0)
        ax.text(-0.12, 1.04, panel, transform=ax.transAxes, fontsize=12, fontweight="bold")
    axes[0].set_ylabel("QMeku interaction energy (kcal/mol)")
    axes[0].legend(frameon=False, fontsize=7, loc="lower right")
    fig.suptitle("Cross-target Kd benchmark: no global correlation", fontsize=11)
    save_figure(fig, "kd_correlation")
    plt.close(fig)


def brd4_figure(data: pd.DataFrame) -> None:
    ids = ["5F63", "3P5O", "4WIV"]
    frame = data.set_index("pdb_id").loc[ids]
    columns = {
        "Experimental": "experimental_dg_kcal_mol",
        "xTB/ALPB": "xtb_e_bind_kcal_mol",
        "MOPAC/COSMO": "mopac_e_bind_kcal_mol",
    }
    colors = {"5F63": "#0072B2", "3P5O": "#D55E00", "4WIV": "#009E73"}
    markers = {"5F63": "o", "3P5O": "s", "4WIV": "^"}
    ranks = pd.DataFrame(index=ids)
    for label, column in columns.items():
        ranks[label] = frame[column].rank(method="min", ascending=True).astype(int)
    fig, ax = plt.subplots(figsize=(5.6, 3.8), layout="constrained")
    xpos = np.arange(len(columns))
    for pdb_id in ids:
        y = ranks.loc[pdb_id].to_numpy()
        ax.plot(xpos, y, color=colors[pdb_id], marker=markers[pdb_id], linewidth=1.8,
                markersize=7, label=pdb_id)
    ax.set_xticks(xpos, list(columns))
    ax.set_yticks([1, 2, 3], ["1 (strongest)", "2", "3 (weakest)"])
    ax.set_ylim(3.25, 0.75)
    ax.set_ylabel("Affinity rank")
    ax.set_title("BRD4 ranking (three crystallographic ligands)")
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    ax.legend(frameon=False, ncols=3, loc="lower center")
    save_figure(fig, "brd4_ranking")
    plt.close(fig)


def write_report(data: pd.DataFrame, metrics: pd.DataFrame) -> None:
    def get(group: str, method: str, kind: str) -> pd.Series:
        return metrics[(metrics.group == group) & (metrics.method == method) & (metrics.correlation == kind)].iloc[0]

    lines = [
        "# QMeku crystal-affinity correlation analysis",
        "",
        "This is an exploratory validation. QMeku values are single-point ligand-pocket interaction energies, not standard binding free energies.",
        "",
        "## Planned strata",
        "",
        "- Primary cross-target panel: six `Kd` complexes in `core_absolute`; one target per complex.",
        "- Expanded `Kd` panel: nine complexes, including repeated BRD4 measurements and the 1BCU stress case.",
        "- BRD4 ranking: three ligands, interpreted descriptively because n=3.",
        "- `Ki` panel: three complexes, interpreted descriptively and not pooled with `Kd`.",
        "- 3PTB remains a development control.",
        "",
        "## Results",
        "",
    ]
    for group, title in (("Kd_core", "Primary cross-target Kd panel"), ("Kd_all", "Expanded Kd panel")):
        lines += [f"### {title}", ""]
        for method in METHODS:
            r = get(group, method, "Pearson")
            rho = get(group, method, "Spearman")
            lines.append(
                f"- {method}: Pearson r={r.estimate:.3f}, 95% bootstrap CI [{r.bootstrap_ci95_low:.3f}, {r.bootstrap_ci95_high:.3f}], "
                f"p={r.p_value_two_sided:.3f}; Spearman ρ={rho.estimate:.3f}, 95% bootstrap CI "
                f"[{rho.bootstrap_ci95_low:.3f}, {rho.bootstrap_ci95_high:.3f}], p={rho.p_value_two_sided:.3f}."
            )
        lines += [""]
    lines += [
        "### BRD4 ranking",
        "",
        "- Experimental order: 5F63 > 3P5O > 4WIV (strongest to weakest).",
        "- xTB reproduces the exact order (Spearman ρ=1.0).",
        "- MOPAC orders 5F63 > 4WIV > 3P5O (Spearman ρ=0.5).",
        "- With n=3 these are descriptive ranking results, not reliable significance tests.",
        "",
        "### Interpretation",
        "",
        "Neither method shows evidence of a general cross-target correlation in this small heterogeneous set. "
        "The xTB result is promising for within-target BRD4 ranking but does not generalize to absolute affinity across targets. "
        "Charged systems and special contexts are influential: 1BCU and 2HB1 receive very favorable xTB interaction energies despite weak experimental affinity, while MOPAC gives positive interaction energies for 2V00, 3COW, and 1AMW.",
        "",
        "## Diagnostics and limitations",
        "",
        "- No observations were excluded after examining the results.",
        "- Pearson and Spearman were both reported because linearity and normal residuals cannot be established credibly with n=6–9.",
        "- Bootstrap intervals use paired resampling with 20,000 draws and seed 20261031; their width reflects the small sample.",
        "- Observations in the expanded panel are not fully independent because BRD4 contributes three complexes.",
        "- Experimental ΔG values derived from Ki are not treated as Kd measurements.",
        "- The benchmark is too small for a stable fitted calibration model.",
        "",
        "## Reproducibility",
        "",
        "Source data: `../results.csv`. Analysis: `../analyze_correlations.py`.",
    ]
    (OUTPUT / "correlation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(INPUT)
    required = ["experimental_dg_kcal_mol", *METHODS.values()]
    if len(data) != 12 or data[required].isna().any().any():
        raise RuntimeError("results.csv must contain twelve complete benchmark rows")
    groups = {
        "Kd_core": (data[(data.affinity_type == "Kd") & (data["set"] == "core_absolute")], True),
        "Kd_all": (data[data.affinity_type == "Kd"], True),
        "BRD4": (data[data.pdb_id.isin(["3P5O", "5F63", "4WIV"])], False),
        "Ki": (data[data.affinity_type == "Ki"], False),
        "All_mixed_descriptive": (data, False),
    }
    rows = []
    for name, (frame, inferential) in groups.items():
        rows.extend(metric_rows(frame, name, inferential))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUTPUT / "correlation_metrics.csv", index=False, float_format="%.8g")
    scatter_figure(data, metrics)
    brd4_figure(data)
    write_report(data, metrics)
    provenance = {
        "source_data": str(INPUT.relative_to(BENCHMARK)),
        "rows": len(data),
        "missing_values": int(data[["experimental_dg_kcal_mol", *METHODS.values()]].isna().sum().sum()),
        "excluded_rows": 0,
        "bootstrap_resamples": N_BOOT,
        "random_seed": SEED,
        "python_packages": {
            "numpy": np.__version__, "pandas": pd.__version__, "scipy": stats.__version__ if hasattr(stats, "__version__") else "see scipy",
            "matplotlib": mpl.__version__,
        },
        "note": "Kd, Ki, and BRD4 strata were specified before calculating correlations in this script.",
    }
    (OUTPUT / "analysis_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
