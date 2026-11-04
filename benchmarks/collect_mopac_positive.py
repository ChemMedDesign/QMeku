#!/usr/bin/env python3
"""Collect benchmark cases with positive MOPAC binding energies for inspection."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path


BENCHMARKS = Path(__file__).resolve().parent
OUTPUT = BENCHMARKS / "mopac_positive_cases"
COPY_FILES = (
    "metadata.json",
    "prepared/receptor.pdb",
    "prepared/ligand.sdf",
    "prepared/preparation.json",
    "prepared/pdb4amber.log",
    "raw/structure.pdb",
    "raw/ligand_crystal.pdb",
    "results/poses.csv",
    "results/qmeku.log",
)


def main() -> None:
    hits: list[tuple[Path, dict[str, str]]] = []
    scanned = 0
    for poses_path in sorted(BENCHMARKS.rglob("poses.csv")):
        if OUTPUT in poses_path.parents:
            continue
        with poses_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                value = row.get("mopac_e_bind", "")
                if not value:
                    continue
                scanned += 1
                if float(value) > 0.0:
                    hits.append((poses_path.parent.parent, row))

    OUTPUT.mkdir(exist_ok=True)
    manifest: list[dict[str, object]] = []
    used_names: set[str] = set()
    for case_dir, row in hits:
        metadata_path = case_dir / "metadata.json"
        metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else {}
        pdb_id = str(metadata.get("pdb_id") or row["pose"].split("_", 1)[0]).upper()
        if pdb_id in used_names:
            raise RuntimeError(f"Duplicate positive case name {pdb_id}: {case_dir}")
        used_names.add(pdb_id)
        destination = OUTPUT / pdb_id
        copied: list[str] = []
        for relative_name in COPY_FILES:
            source = case_dir / relative_name
            if not source.is_file():
                continue
            target = destination / relative_name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            copied.append(relative_name)

        manifest.append({
            "pdb_id": pdb_id,
            "pose": row.get("pose", ""),
            "source_case": case_dir.relative_to(BENCHMARKS.parent).as_posix(),
            "ligand_ccd": metadata.get("ligand_ccd", ""),
            "ligand_charge": row.get("ligand_charge", ""),
            "pocket_charge": row.get("pocket_charge", ""),
            "pocket_residues": row.get("pocket_residues", ""),
            "pocket_atoms": row.get("pocket_atoms", ""),
            "solvent": row.get("solvent", ""),
            "xtb_e_bind_kcal_mol": row.get("xtb_e_bind", ""),
            "mopac_e_bind_kcal_mol": row.get("mopac_e_bind", ""),
            "copied_files": ";".join(copied),
        })

    manifest.sort(key=lambda item: float(item["mopac_e_bind_kcal_mol"]), reverse=True)
    manifest_path = OUTPUT / "manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    print(f"Scanned {scanned} MOPAC results")
    print(f"Collected {len(manifest)} positive cases in {OUTPUT}")
    for item in manifest:
        print(f"  {item['pdb_id']}: {float(item['mopac_e_bind_kcal_mol']):+.3f} kcal/mol")


if __name__ == "__main__":
    main()
