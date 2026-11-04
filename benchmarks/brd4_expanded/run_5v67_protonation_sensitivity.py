#!/usr/bin/env python3
"""Recalculate 5V67 with the alternate protonated piperazine nitrogen."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent
ORIGINAL = BENCHMARK / "kd_itc_series" / "5V67"
OUTPUT = BENCHMARK / "protonation_sensitivity" / "5V67_alternate_piperazine_N"
SHARED = BENCHMARK.parent / "affinity_crystal" / "prepare_benchmark.py"


def load_shared():
    spec = importlib.util.spec_from_file_location("affinity_crystal_preparation", SHARED)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SHARED}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def heavy_coordinates(molecule: Chem.Mol):
    conformer = molecule.GetConformer()
    return [tuple(conformer.GetAtomPosition(atom.GetIdx())) for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]


def one_energy(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle))


def main() -> int:
    source = Chem.MolFromMolFile(str(ORIGINAL / "raw" / "ligand_crystal_unprepared.sdf"), removeHs=False)
    if source is None:
        raise RuntimeError("could not read 5V67 crystal ligand")
    old_heavy = heavy_coordinates(source)
    molecule = Chem.RemoveHs(source)
    # In the deposited atom order, atoms 39 and 40 are the two piperazine
    # nitrogens. The primary preparation protonates 39; this sensitivity case
    # protonates the opposite site, atom 40.
    alternate = molecule.GetAtomWithIdx(40)
    if alternate.GetAtomicNum() != 7 or alternate.GetDegree() != 3 or alternate.GetFormalCharge() != 0:
        raise RuntimeError("atom 40 is no longer the expected neutral tertiary piperazine nitrogen")
    alternate.SetFormalCharge(1)
    alternate.SetNumExplicitHs(1)
    alternate.SetNoImplicit(True)
    Chem.SanitizeMol(molecule)
    prepared = Chem.AddHs(molecule, addCoords=True)
    prepared.GetConformer().Set3D(True)
    prepared.SetProp("_Name", "5V67_IBI_crystal_pH7.4_alternate_piperazine_N40")
    new_heavy = heavy_coordinates(prepared)
    rmsd = math.sqrt(
        sum(sum((a - b) ** 2 for a, b in zip(left, right)) for left, right in zip(old_heavy, new_heavy))
        / len(new_heavy)
    )
    if rmsd > 1e-6 or Chem.GetFormalCharge(prepared) != 1:
        raise RuntimeError(f"invalid alternate protomer: charge={Chem.GetFormalCharge(prepared)}, RMSD={rmsd}")
    ligand = OUTPUT / "ligand.sdf"
    ligand.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(ligand))
    writer.write(prepared)
    writer.close()
    receptor = ORIGINAL / "prepared" / "receptor.pdb"
    validation = load_shared().validate_for_qmeku(receptor, ligand)
    outdir = OUTPUT / "results"
    command = [
        str(ROOT / "bin" / "qmeku"), str(receptor), str(ligand),
        "--outdir", str(outdir), "--method", "both", "--total-only",
        "--no-gas", "--workers", "2", "--quiet",
    ]
    completed = subprocess.run(command, cwd=ROOT, check=False)
    if completed.returncode:
        raise RuntimeError(f"QMeku failed ({completed.returncode})")
    old = one_energy(ORIGINAL / "results" / "poses.csv")
    new = one_energy(outdir / "poses.csv")
    comparison = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "primary_protonated_atom": 39,
        "alternate_protonated_atom": 40,
        "ligand_charge": 1,
        "heavy_atom_rmsd_A": rmsd,
        "qmeku_validation": validation,
        "primary_xtb": float(old["xtb_e_bind"]),
        "alternate_xtb": float(new["xtb_e_bind"]),
        "delta_xtb": float(new["xtb_e_bind"]) - float(old["xtb_e_bind"]),
        "primary_mopac": float(old["mopac_e_bind"]),
        "alternate_mopac": float(new["mopac_e_bind"]),
        "delta_mopac": float(new["mopac_e_bind"]) - float(old["mopac_e_bind"]),
        "qmeku_command": command,
    }
    (OUTPUT / "comparison.json").write_text(json.dumps(comparison, indent=2, sort_keys=True) + "\n")
    print(OUTPUT / "comparison.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
