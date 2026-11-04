#!/usr/bin/env python3
"""Collect one QMeku energy row per expanded BRD4 complex."""

from __future__ import annotations

import csv
from pathlib import Path


BENCHMARK = Path(__file__).resolve().parent
ROOT = BENCHMARK.parents[1]


with (BENCHMARK / "manifest.csv").open(newline="", encoding="utf-8") as handle:
    manifest = list(csv.DictReader(handle))

fieldnames = [
    "pdb_id", "set", "ligand_name", "ligand_ccd", "activity_type", "activity_value",
    "activity_error", "activity_unit", "p_activity", "experimental_dg_kcal_mol",
    "xtb_e_bind_kcal_mol", "mopac_e_bind_kcal_mol", "pocket_residues", "pocket_atoms",
    "pocket_charge", "ligand_charge", "status",
]
rows: list[dict[str, str]] = []
for item in manifest:
    result = BENCHMARK / item["set"] / item["pdb_id"] / "results" / "poses.csv"
    receptor = ROOT / item["receptor_path"]
    ligand = ROOT / item["ligand_path"]
    output = {key: item.get(key, "") for key in fieldnames}
    stale = result.is_file() and any(
        path.is_file() and path.stat().st_mtime > result.stat().st_mtime
        for path in (receptor, ligand)
    )
    if not result.is_file() or stale:
        output["status"] = "pending"
        rows.append(output)
        continue
    with result.open(newline="", encoding="utf-8") as handle:
        energy = next(csv.DictReader(handle))
    output.update({
        "xtb_e_bind_kcal_mol": energy.get("xtb_e_bind", ""),
        "mopac_e_bind_kcal_mol": energy.get("mopac_e_bind", ""),
        "pocket_residues": energy.get("pocket_residues", ""),
        "pocket_atoms": energy.get("pocket_atoms", ""),
        "pocket_charge": energy.get("pocket_charge", ""),
        "ligand_charge": energy.get("ligand_charge", ""),
        "status": "completed",
    })
    item["status"] = "completed"
    rows.append(output)

with (BENCHMARK / "results.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

temporary = BENCHMARK / "manifest.csv.part"
with temporary.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
    writer.writeheader()
    writer.writerows(manifest)
temporary.replace(BENCHMARK / "manifest.csv")
