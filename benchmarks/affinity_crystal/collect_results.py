#!/usr/bin/env python3
"""Collect one total-energy row per completed benchmark complex."""

from __future__ import annotations

import csv
from pathlib import Path


BENCHMARK = Path(__file__).resolve().parent
ROOT = BENCHMARK.parents[1]


with (BENCHMARK / "manifest.csv").open(newline="", encoding="utf-8") as handle:
    manifest = list(csv.DictReader(handle))

fieldnames = [
    "pdb_id", "set", "ligand_ccd", "affinity_type", "affinity_value",
    "affinity_unit", "experimental_dg_kcal_mol", "xtb_e_bind_kcal_mol",
    "mopac_e_bind_kcal_mol", "pocket_residues", "pocket_atoms",
    "pocket_charge", "ligand_charge", "status",
]
rows = []
for item in manifest:
    result = BENCHMARK / item["set"] / item["pdb_id"] / "results" / "poses.csv"
    receptor = ROOT / item["receptor_path"]
    ligand = ROOT / item["ligand_path"]
    output = {key: item.get(key, "") for key in fieldnames}
    stale = result.is_file() and any(
        p.is_file() and p.stat().st_mtime > result.stat().st_mtime
        for p in (receptor, ligand)
    )
    if not result.is_file() or stale:
        output["status"] = "pending"
        if stale and receptor.is_file() and ligand.is_file():
            item["status"] = "prepared"
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

if manifest:
    tmp_manifest = BENCHMARK / "manifest.csv.part"
    with tmp_manifest.open("w", newline="", encoding="utf-8") as handle:
        mwriter = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        mwriter.writeheader()
        mwriter.writerows(manifest)
    tmp_manifest.replace(BENCHMARK / "manifest.csv")
