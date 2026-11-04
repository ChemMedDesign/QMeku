#!/usr/bin/env python3
"""Prepare BRD4 receptors and crystallographic ligands for QMeku.

Protein preparation delegates to the audited affinity-crystal workflow. Ligand
heavy atoms and their coordinates are preserved; only hydrogens and the explicit
benchmark protonation state at pH 7.4 are added.
"""

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
        raise RuntimeError(f"could not load shared preparation module: {SHARED}")
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


def is_carbonyl_carbon(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() == 6 and any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetOtherAtom(atom).GetAtomicNum() in (8, 16)
        for bond in atom.GetBonds()
    )


def protonate_basic_tertiary_amine(molecule: Chem.Mol) -> tuple[int, list[dict[str, object]]]:
    candidates: list[tuple[tuple[int, int, int, int], Chem.Atom]] = []
    audit: list[dict[str, object]] = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 0 or atom.GetIsAromatic() or atom.GetDegree() != 3:
            continue
        neighbors = list(atom.GetNeighbors())
        if any(is_carbonyl_carbon(neighbor) for neighbor in neighbors):
            continue
        methyl_neighbors = sum(
            neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1
            for neighbor in neighbors
        )
        aromatic_neighbors = sum(neighbor.GetIsAromatic() for neighbor in neighbors)
        hetero_neighbors = sum(neighbor.GetAtomicNum() not in (1, 6) for neighbor in neighbors)
        carbon_neighbors = sum(neighbor.GetAtomicNum() == 6 for neighbor in neighbors)
        # Prefer the N-methyl terminus of substituted piperazines. Avoid anilines
        # and other conjugated nitrogens when a more aliphatic site is available.
        score = (methyl_neighbors, -aromatic_neighbors, carbon_neighbors, -atom.GetIdx())
        candidates.append((score, atom))
        audit.append({
            "atom_index": atom.GetIdx(),
            "methyl_neighbors": methyl_neighbors,
            "aromatic_neighbors": aromatic_neighbors,
            "hetero_neighbors": hetero_neighbors,
            "score": list(score),
        })
    if not candidates:
        raise RuntimeError("no neutral tertiary aliphatic amine available for the requested +1 protomer")
    candidates.sort(key=lambda item: item[0], reverse=True)
    atom = candidates[0][1]
    atom.SetFormalCharge(1)
    atom.SetNumExplicitHs(1)
    atom.SetNoImplicit(True)
    return atom.GetIdx(), audit


def prepare_ligand(case: Path, row: dict[str, str], force: bool) -> tuple[Path, dict[str, object]]:
    pdb_id = row["pdb_id"].upper()
    ccd = row["ligand_ccd"].upper()
    expected_charge = int(row["expected_charge"])
    output = case / "prepared" / "ligand.sdf"
    if output.is_file() and not force:
        return output, {"reused": True}
    output.parent.mkdir(parents=True, exist_ok=True)
    source = case / "raw" / "ligand_crystal_unprepared.sdf"
    original = Chem.MolFromMolFile(str(source), removeHs=False)
    if original is None or original.GetNumConformers() != 1:
        raise RuntimeError("RDKit could not read one crystallographic ligand conformer")
    original_heavy = heavy_coordinates(original)
    molecule = Chem.RemoveHs(original)
    initial_charge = Chem.GetFormalCharge(molecule)
    changes: list[str] = []
    candidate_audit: list[dict[str, object]] = []
    if initial_charge != expected_charge:
        if row["protonation_rule"] != "basic_tertiary_amine" or expected_charge != initial_charge + 1:
            raise RuntimeError(
                f"CCD charge {initial_charge:+d} cannot be changed to expected {expected_charge:+d} "
                f"with rule {row['protonation_rule']!r}"
            )
        atom_index, candidate_audit = protonate_basic_tertiary_amine(molecule)
        changes.append(f"protonated tertiary aliphatic N atom {atom_index}")
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistryFrom3D(molecule)
    prepared = Chem.AddHs(molecule, addCoords=True)
    prepared.GetConformer().Set3D(True)
    final_charge = Chem.GetFormalCharge(prepared)
    if final_charge != expected_charge:
        raise RuntimeError(f"prepared charge {final_charge:+d}; expected {expected_charge:+d}")
    new_heavy = heavy_coordinates(prepared)
    if len(original_heavy) != len(new_heavy):
        raise RuntimeError("ligand heavy atom count changed during preparation")
    rmsd = math.sqrt(
        sum(sum((a - b) ** 2 for a, b in zip(left, right)) for left, right in zip(original_heavy, new_heavy))
        / len(new_heavy)
    )
    if rmsd > 1e-6:
        raise RuntimeError(f"ligand heavy coordinates changed: RMSD={rmsd}")
    prepared.SetProp("_Name", f"{pdb_id}_{ccd}_crystal_pH{PH}")
    prepared.SetProp("PDB_ID", pdb_id)
    prepared.SetProp("CCD_ID", ccd)
    prepared.SetProp("PROTONATION_PH", str(PH))
    prepared.SetProp("PROTONATION_NOTES", "; ".join(changes) if changes else "CCD protomer retained")
    writer = Chem.SDWriter(str(output))
    writer.write(prepared)
    writer.close()
    return output, {
        "reused": False,
        "initial_ccd_charge": initial_charge,
        "charge": final_charge,
        "protonation_rule": row["protonation_rule"],
        "protonation_changes": changes or ["CCD protomer retained"],
        "protonation_candidate_audit": candidate_audit,
        "atom_count": prepared.GetNumAtoms(),
        "heavy_atom_count": sum(atom.GetAtomicNum() > 1 for atom in prepared.GetAtoms()),
        "heavy_atom_rmsd_A": rmsd,
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
            (case / "prepared" / "preparation.json").write_text(
                json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            row["status"] = "prepared"
            completed.add(pdb_id)
            print(
                f"  ready: ligand {validation['ligand_charge']:+d}; "
                f"pocket {validation['pocket_residues_at_4A']} residues / {validation['pocket_atoms']} atoms",
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
