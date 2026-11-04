#!/usr/bin/env python3
"""Reconcile the neutral BRD4 manifest with completed QMEku calculations."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.csv"
OUTPUT = ROOT / "results.csv"


def read_single_pose(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one pose in {path}, found {len(rows)}")
    return rows[0]


def heavy_atom_count(sdf_path: Path) -> int:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    mol = next((item for item in supplier if item is not None), None)
    if mol is None:
        raise ValueError(f"RDKit could not read {sdf_path}")
    return sum(atom.GetAtomicNum() > 1 for atom in mol.GetAtoms())


def main() -> None:
    with MANIFEST.open(newline="") as handle:
        manifest = list(csv.DictReader(handle))
        manifest_fields = list(manifest[0])

    output_rows: list[dict[str, object]] = []
    failures: list[str] = []
    for row in manifest:
        pdb_id = row["pdb_id"]
        case_dir = ROOT / row["set"] / pdb_id
        pose_path = case_dir / "results" / "poses.csv"
        prep_path = case_dir / "prepared" / "preparation.json"
        try:
            pose = read_single_pose(pose_path)
            prep = json.loads(prep_path.read_text())
            ic50_nm = float(row["ic50_geomean_nM"])
            exact_pic50 = 9.0 - math.log10(ic50_nm)
            ligand = prep["ligand"]
            validation = prep["qmeku_validation"]
            result = {
                "pdb_id": pdb_id,
                "set": row["set"],
                "ligand_ccd": row["ligand_ccd"],
                "ligand_name": row["ligand_name"],
                "assay_series": row["assay_series"],
                "ic50_geomean_nM": f"{ic50_nm:.8g}",
                "pIC50_exact": f"{exact_pic50:.8f}",
                "resolution_A": row["resolution_A"],
                "ligand_charge": pose["ligand_charge"],
                "ligand_heavy_atoms": ligand.get(
                    "heavy_atom_count", heavy_atom_count(case_dir / "prepared" / "ligand.sdf")
                ),
                "pocket_residues": pose["pocket_residues"],
                "pocket_atoms": pose["pocket_atoms"],
                "pocket_charge": pose["pocket_charge"],
                "even_electron_model": validation["even_electron_model"],
                "xtb_e_bind_kcal_mol": pose["xtb_e_bind"],
                "mopac_e_bind_kcal_mol": pose["mopac_e_bind"],
                "qmeku_seconds": pose["seconds"],
                "source_doi": row["source_doi"],
            }
            for method in ("xtb", "mopac"):
                if not result[f"{method}_e_bind_kcal_mol"]:
                    raise ValueError(f"missing {method} binding energy")
            output_rows.append(result)
            row["status"] = "completed_both"
        except Exception as exc:  # retain all failures for one actionable report
            failures.append(f"{pdb_id}: {exc}")
            row["status"] = "result_error"

    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields)
        writer.writeheader()
        writer.writerows(manifest)

    if failures:
        raise RuntimeError("Result reconciliation failed:\n" + "\n".join(failures))

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} complete cases to {OUTPUT}")


if __name__ == "__main__":
    main()
