#!/usr/bin/env python3
"""Calculate panel-specific correlations for the expanded BRD4 benchmark."""

from __future__ import annotations

import csv
from pathlib import Path

from scipy.stats import pearsonr, spearmanr


BENCHMARK = Path(__file__).resolve().parent


with (BENCHMARK / "results.csv").open(newline="", encoding="utf-8") as handle:
    results = [row for row in csv.DictReader(handle) if row["status"] == "completed"]

analyses = [
    ("ic50_primary", [row for row in results if row["set"] == "ic50_series"], "none"),
    ("kd_itc_primary", [row for row in results if row["set"] == "kd_itc_series"], "none"),
    (
        "kd_itc_sensitivity_without_3MXF",
        [row for row in results if row["set"] == "kd_itc_series" and row["pdb_id"] != "3MXF"],
        "3MXF excluded because xTB/MOPAC disagree in sign; sensitivity analysis only",
    ),
]

fieldnames = [
    "analysis", "method", "experimental_variable", "expected_sign", "n", "excluded",
    "pearson_r", "pearson_p", "spearman_rho", "spearman_p",
]
metrics: list[dict[str, object]] = []
for analysis, rows, excluded in analyses:
    for method in ("xtb", "mopac"):
        calculated = [float(row[f"{method}_e_bind_kcal_mol"]) for row in rows]
        experimental = [float(row["p_activity"]) for row in rows]
        pearson = pearsonr(calculated, experimental)
        spearman = spearmanr(calculated, experimental)
        metrics.append({
            "analysis": analysis,
            "method": method,
            "experimental_variable": "pIC50" if analysis.startswith("ic50") else "pKd",
            "expected_sign": "negative",
            "n": len(rows),
            "excluded": excluded,
            "pearson_r": f"{pearson.statistic:.6f}",
            "pearson_p": f"{pearson.pvalue:.6f}",
            "spearman_rho": f"{spearman.statistic:.6f}",
            "spearman_p": f"{spearman.pvalue:.6f}",
        })

with (BENCHMARK / "correlation_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(metrics)

print(f"Wrote {len(metrics)} panel/method analyses to correlation_metrics.csv")
