#!/usr/bin/env python3
"""Prepare neutral BRD4 crystal complexes for QMeku."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent
MANIFEST = BENCHMARK / "manifest.csv"
SHARED = BENCHMARK.parent / "affinity_crystal" / "prepare_benchmark.py"
PH = 7.4


def load_shared():
    spec = importlib.util.spec_from_file_location("affinity_crystal_preparation", SHARED)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SHARED}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.BENCHMARK = BENCHMARK
    module.MANIFEST = MANIFEST
    return module


def heavy_coordinates(molecule: Chem.Mol) -> list[tuple[float, float, float]]:
    conformer = molecule.GetConformer()
    return [
        tuple(conformer.GetAtomPosition(atom.GetIdx()))
        for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1
    ]


def prepare_ligand(case: Path, row: dict[str, str], force: bool) -> tuple[Path, dict[str, object]]:
    pdb_id = row["pdb_id"].upper()
    ccd = row["ligand_ccd"].upper()
    expected_charge = int(row["expected_charge"])
    output = case / "prepared" / "ligand.sdf"
    if output.is_file() and not force:
        return output, {"reused": True}
    output.parent.mkdir(parents=True, exist_ok=True)
    source = Chem.MolFromMolFile(str(case / "raw" / "ligand_crystal_unprepared.sdf"), removeHs=False)
    if source is None or source.GetNumConformers() != 1:
        raise RuntimeError("RDKit could not read one crystallographic ligand conformer")
    original_heavy = heavy_coordinates(source)
    molecule = Chem.RemoveHs(source)
    initial_charge = Chem.GetFormalCharge(molecule)
    if initial_charge != expected_charge:
        raise RuntimeError(f"CCD-derived ligand charge {initial_charge:+d}; expected {expected_charge:+d}")
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistryFrom3D(molecule)
    prepared = Chem.AddHs(molecule, addCoords=True)
    prepared.GetConformer().Set3D(True)
    final_charge = Chem.GetFormalCharge(prepared)
    if final_charge != expected_charge:
        raise RuntimeError(f"prepared ligand charge {final_charge:+d}; expected {expected_charge:+d}")
    new_heavy = heavy_coordinates(prepared)
    if len(original_heavy) != len(new_heavy):
        raise RuntimeError("heavy atom count changed during ligand preparation")
    rmsd = math.sqrt(
        sum(sum((a - b) ** 2 for a, b in zip(left, right)) for left, right in zip(original_heavy, new_heavy))
        / len(new_heavy)
    )
    if rmsd > 1e-6:
        raise RuntimeError(f"ligand heavy coordinates changed: RMSD={rmsd}")
    prepared.SetProp("_Name", f"{pdb_id}_{ccd}_crystal_neutral_pH{PH}")
    prepared.SetProp("PDB_ID", pdb_id)
    prepared.SetProp("CCD_ID", ccd)
    prepared.SetProp("PROTONATION_PH", str(PH))
    prepared.SetProp("PROTONATION_NOTES", "neutral CCD-derived protomer retained")
    writer = Chem.SDWriter(str(output))
    writer.write(prepared)
    writer.close()
    return output, {
        "reused": False,
        "initial_ccd_charge": initial_charge,
        "charge": final_charge,
        "atom_count": prepared.GetNumAtoms(),
        "heavy_atom_count": len(new_heavy),
        "heavy_atom_rmsd_A": rmsd,
        "protonation": "neutral CCD-derived protomer retained",
    }


def write_manifest(rows: list[dict[str, str]]) -> None:
    temporary = MANIFEST.with_suffix(".csv.part")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(MANIFEST)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", help="optional PDB IDs")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    wanted = {item.upper() for item in args.ids} if args.ids else None
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = [row for row in rows if wanted is None or row["pdb_id"].upper() in wanted]
    shared = load_shared()
    amber_bin = shared.ambertools_bin()
    completed: set[str] = set()
    failures: list[tuple[str, str]] = []
    for index, row in enumerate(selected, start=1):
        pdb_id = row["pdb_id"].upper()
        case = BENCHMARK / row["set"] / pdb_id
        print(f"[{index}/{len(selected)}] {pdb_id}", file=sys.stderr)
        try:
            receptor, receptor_meta = shared.prepare_receptor(case, pdb_id, amber_bin, args.force)
            ligand, ligand_meta = prepare_ligand(case, row, args.force)
            validation = shared.validate_for_qmeku(receptor, ligand)
            preparation_path = case / "prepared" / "preparation.json"
            if receptor_meta.get("reused") and ligand_meta.get("reused") and preparation_path.is_file():
                metadata = json.loads(preparation_path.read_text(encoding="utf-8"))
                metadata["revalidated_utc"] = datetime.now(timezone.utc).isoformat()
                metadata["qmeku_validation"] = validation
            else:
                metadata = {
                    "pdb_id": pdb_id,
                    "prepared_utc": datetime.now(timezone.utc).isoformat(),
                    "protein_tool": "pdb4amber --reduce (AmberTools)",
                    "ligand_tool": f"RDKit {Chem.rdBase.rdkitVersion}",
                    "protonation_ph": PH,
                    "receptor": receptor_meta,
                    "ligand": ligand_meta,
                    "qmeku_validation": validation,
                }
            preparation_path.write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            row["status"] = "prepared"
            completed.add(pdb_id)
            print(
                f"  ready: q={validation['ligand_charge']:+d}; "
                f"pocket={validation['pocket_residues_at_4A']} residues/{validation['pocket_atoms']} atoms",
                file=sys.stderr,
            )
        except Exception as exc:
            failures.append((pdb_id, str(exc)))
            print(f"  FAILED: {exc}", file=sys.stderr)
    write_manifest(rows)
    print(
        f"Preparation reconciliation: requested={len(selected)}, prepared={len(completed)}, failed={len(failures)}",
        file=sys.stderr,
    )
    for pdb_id, message in failures:
        print(f"  {pdb_id}: {message}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
