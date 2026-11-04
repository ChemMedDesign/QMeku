#!/usr/bin/env python3
"""Prepare downloaded crystal receptors and ligands for QMeku.

Protein preparation uses AmberTools ``pdb4amber --reduce`` at its default
protonation model. Ligands retain every crystallographic heavy-atom coordinate;
ionizable groups are assigned to the documented pH 7.4 protomer and hydrogens are
added with RDKit. No heavy atom is optimized.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent
MANIFEST = BENCHMARK / "manifest.csv"
PH = 7.4
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Explicit benchmark protomer choices at pH 7.4. These are not inferred later
# from the total charge: the transformations below identify the chemical group.
LIGAND_CHARGES = {
    "2V00": 0,   # V15
    "1ADL": -1,  # arachidonate
    "3P5O": 0,   # I-BET
    "3COW": 0,   # 52H
    "1AMW": -3,  # ADP3-
    "1BM7": -1,  # flufenamate
    "5F63": 1,   # protonated N-methylpiperidine
    "4WIV": 0,   # UMB32
    "1E66": 1,   # huprine X, protonated ring nitrogen
    "2HB1": -2,  # dicarboxylate
    "1BCU": 1,   # proflavine, protonated acridine nitrogen
}

CARBOXYLATE_COUNTS = {"1ADL": 1, "1BM7": 1, "2HB1": 2}
AROMATIC_N_PROTONATION = {"1E66", "1BCU"}
ALIPHATIC_N_PROTONATION = {"5F63"}


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def ambertools_bin() -> Path:
    candidates: list[Path] = []
    if os.environ.get("AMBERTOOLS_BIN"):
        candidates.append(Path(os.environ["AMBERTOOLS_BIN"]))
    candidates.append(Path(sys.executable).resolve().parent)
    executable = shutil.which("pdb4amber")
    if executable:
        candidates.append(Path(executable).resolve().parent)
    candidates.extend([
        Path.home() / "miniforge3" / "envs" / "mddeku-md" / "bin",
        Path.home() / "miniforge3" / "envs" / "mdeku" / "bin",
    ])
    for directory in candidates:
        try:
            p4a = directory / "pdb4amber"
            red = directory / "reduce"
            if p4a.is_file() and os.access(p4a, os.X_OK) and red.is_file() and os.access(red, os.X_OK):
                return directory
        except OSError:
            continue
    raise RuntimeError("AmberTools with pdb4amber and reduce was not found (put it in PATH or set AMBERTOOLS_BIN)")


def selected_chains(case: Path) -> set[str]:
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    return {" " if chain == "_" else chain for chain in metadata["selected_protein_chains_within_8A"]}


def mark_geometric_chain_breaks(path: Path) -> list[str]:
    """Insert TER where consecutive residue records are not peptide bonded and renumber sequentially
    without insertion codes so pdb4amber/reduce never misidentifies insertion codes as chain termini."""
    atom_lines = [line for line in path.read_text().splitlines() if line.startswith("ATOM  ")]
    residues: list[tuple[tuple[str, str], list[str]]] = []
    for line in atom_lines:
        key = (line[21:22], line[22:27])
        if not residues or residues[-1][0] != key:
            residues.append((key, []))
        residues[-1][1].append(line)

    output: list[str] = []
    breaks: list[str] = []
    previous_key = None
    previous_atoms: list[str] | None = None
    for seq_num, (key, atoms) in enumerate(residues, start=1):
        if previous_atoms is not None:
            previous_c = next((line for line in previous_atoms if line[12:16].strip() == "C"), None)
            current_n = next((line for line in atoms if line[12:16].strip() == "N"), None)
            distance = None
            if previous_c is not None and current_n is not None and previous_key[0] == key[0]:
                a = tuple(float(previous_c[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
                b = tuple(float(current_n[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
                distance = math.dist(a, b)
            if previous_key[0] != key[0] or distance is None or distance >= 2.0:
                output.append("TER")
                label = f"{previous_key[0] or '_'}:{previous_key[1].strip()} -> {key[0] or '_'}:{key[1].strip()}"
                breaks.append(label + (f" ({distance:.2f} A)" if distance is not None else ""))
        output.extend(f"{line[:22]}{seq_num:4d} {line[27:]}" for line in atoms)
        previous_key, previous_atoms = key, atoms
    output.extend(["TER", "END"])
    path.write_text("\n".join(output) + "\n", encoding="ascii")
    return breaks


def protein_input(case: Path, pdb_id: str, destination: Path) -> list[str]:
    """Create the heavy-atom protein passed to pdb4amber.

    1ADL needs its polymeric MHO40 and CSD117 restored from HETATM records and
    reverted to their standard parents. Extra oxidation oxygens are removed.
    """
    if pdb_id != "1ADL":
        shutil.copyfile(case / "raw" / "receptor_protein_only_unprepared.pdb", destination)
        return ["geometric chain break: " + value for value in mark_geometric_chain_breaks(destination)]

    chains = selected_chains(case)
    source = (case / "raw" / "structure.pdb").read_text(encoding="ascii", errors="replace").splitlines()
    allowed = {
        ("A", "40", "MHO"): ("MET", {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"}),
        ("A", "117", "CSD"): ("CYS", {"N", "CA", "C", "O", "CB", "SG"}),
    }
    output: list[str] = []
    for line in source:
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or line[21:22] not in chains or line[16:17] not in (" ", "A"):
            continue
        key = (line[21:22], line[22:26].strip(), line[17:20].strip())
        if line.startswith("ATOM  "):
            output.append(line[:16] + " " + line[17:])
        elif key in allowed:
            parent, atom_names = allowed[key]
            if line[12:16].strip() in atom_names:
                output.append("ATOM  " + line[6:17] + f"{parent:>3}" + line[20:])
    output.append("END")
    destination.write_text("\n".join(output) + "\n", encoding="ascii")
    changes = ["MHO A40 reverted to MET; sulfoxide OD1 removed", "CSD A117 reverted to CYS; OD1/OD2 removed"]
    changes.extend("geometric chain break: " + value for value in mark_geometric_chain_breaks(destination))
    return changes


def prepare_receptor(case: Path, pdb_id: str, amber_bin: Path, force: bool) -> tuple[Path, dict[str, object]]:
    output = case / "prepared" / "receptor.pdb"
    log_output = case / "prepared" / "pdb4amber.log"
    if output.is_file() and not force:
        return output, {"reused": True}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"qmeku_{pdb_id}_") as tmp_name:
        tmp = Path(tmp_name)
        heavy = tmp / "protein_heavy.pdb"
        modifications = protein_input(case, pdb_id, heavy)
        amber_output = tmp / "receptor_amber.pdb"
        amber_log = tmp / "pdb4amber.log"
        env = os.environ.copy()
        env["PATH"] = str(amber_bin) + os.pathsep + env.get("PATH", "")
        command = [
            str(amber_bin / "pdb4amber"), "-i", str(heavy), "-o", str(amber_output),
            "-p", "-a", "-d", "--reduce", "--add-missing-atoms", "-l", str(amber_log),
        ]
        completed = subprocess.run(command, cwd=tmp, env=env, text=True, capture_output=True, check=False)
        if not amber_output.is_file() or amber_output.stat().st_size == 0:
            raise RuntimeError(f"pdb4amber failed ({completed.returncode}): {completed.stderr[-500:]}")
        atoms = [line for line in amber_output.read_text().splitlines() if line.startswith("ATOM  ")]
        output.write_text("\n".join(atoms) + "\nEND\n", encoding="ascii")
        log_text = amber_log.read_text(errors="replace") if amber_log.is_file() else ""
        log_output.write_text(log_text, encoding="utf-8")

    residue_names = {line[17:20].strip() for line in atoms}
    generic_his = sorted({(line[21:22], line[22:26].strip()) for line in atoms if line[17:20].strip() == "HIS"})
    hydrogen_count = sum(line[76:78].strip().upper() == "H" for line in atoms)
    if not atoms or hydrogen_count == 0 or generic_his:
        raise RuntimeError(f"invalid prepared receptor: atoms={len(atoms)}, H={hydrogen_count}, HIS={generic_his}")
    return output, {
        "reused": False,
        "pdb4amber_command": command,
        "pdb4amber_returncode": completed.returncode,
        "heavy_atom_modifications": modifications,
        "atom_count": len(atoms),
        "hydrogen_count": hydrogen_count,
        "residue_names": sorted(residue_names),
    }


def set_carboxylates(mol: Chem.Mol, expected: int) -> list[int]:
    changed: list[int] = []
    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6:
            continue
        oxygens = []
        for bond in carbon.GetBonds():
            other = bond.GetOtherAtom(carbon)
            if other.GetAtomicNum() == 8:
                oxygens.append((bond, other))
        if len(oxygens) != 2 or not any(b.GetBondType() == Chem.BondType.DOUBLE for b, _ in oxygens):
            continue
        singles = [atom for bond, atom in oxygens if bond.GetBondType() == Chem.BondType.SINGLE]
        if len(singles) == 1:
            oxygen = singles[0]
            oxygen.SetFormalCharge(-1)
            oxygen.SetNumExplicitHs(0)
            oxygen.SetNoImplicit(True)
            changed.append(oxygen.GetIdx())
    if len(changed) != expected:
        raise RuntimeError(f"expected {expected} carboxylates, found {len(changed)}")
    return changed


def set_adp_phosphates(mol: Chem.Mol) -> list[int]:
    changed: list[int] = []
    for phosphorus in (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15):
        for bond in phosphorus.GetBonds():
            oxygen = bond.GetOtherAtom(phosphorus)
            if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 1 or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            oxygen.SetFormalCharge(-1)
            oxygen.SetNumExplicitHs(0)
            oxygen.SetNoImplicit(True)
            changed.append(oxygen.GetIdx())
    if len(changed) != 3:
        raise RuntimeError(f"expected three deprotonated terminal ADP oxygens, found {len(changed)}")
    return changed


def protonate_aromatic_n(mol: Chem.Mol) -> int:
    candidates = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetIsAromatic() and atom.GetDegree() == 2]
    if len(candidates) != 1:
        raise RuntimeError(f"expected one protonatable aromatic N, found {[a.GetIdx() for a in candidates]}")
    atom = candidates[0]
    atom.SetFormalCharge(1)
    atom.SetNumExplicitHs(1)
    atom.SetNoImplicit(True)
    return atom.GetIdx()


def protonate_aliphatic_n(mol: Chem.Mol) -> int:
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetIsAromatic() or atom.GetDegree() != 3:
            continue
        is_amide = any(
            neighbor.GetAtomicNum() == 6
            and any(bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(neighbor).GetAtomicNum() == 8
                    for bond in neighbor.GetBonds())
            for neighbor in atom.GetNeighbors()
        )
        if not is_amide and atom.GetTotalNumHs() == 0:
            candidates.append(atom)
    if len(candidates) != 1:
        raise RuntimeError(f"expected one tertiary basic N, found {[a.GetIdx() for a in candidates]}")
    atom = candidates[0]
    atom.SetFormalCharge(1)
    atom.SetNumExplicitHs(1)
    atom.SetNoImplicit(True)
    return atom.GetIdx()


def heavy_coordinates(mol: Chem.Mol) -> list[tuple[float, float, float]]:
    conf = mol.GetConformer()
    return [tuple(conf.GetAtomPosition(atom.GetIdx())) for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]


def prepare_ligand(case: Path, pdb_id: str, ccd: str, force: bool) -> tuple[Path, dict[str, object]]:
    output = case / "prepared" / "ligand.sdf"
    if output.is_file() and not force:
        return output, {"reused": True}
    source = case / "raw" / "ligand_crystal_unprepared.sdf"
    original = Chem.MolFromMolFile(str(source), removeHs=False)
    if original is None or original.GetNumConformers() != 1:
        raise RuntimeError("RDKit could not read one crystallographic ligand conformer")
    original_heavy = heavy_coordinates(original)
    mol = Chem.RemoveHs(original)
    changes: list[str] = []
    if pdb_id in CARBOXYLATE_COUNTS:
        atoms = set_carboxylates(mol, CARBOXYLATE_COUNTS[pdb_id])
        changes.append(f"deprotonated carboxylate O atoms {atoms}")
    if pdb_id == "1AMW":
        atoms = set_adp_phosphates(mol)
        changes.append(f"ADP terminal phosphate O atoms deprotonated {atoms}")
    if pdb_id in AROMATIC_N_PROTONATION:
        atom = protonate_aromatic_n(mol)
        changes.append(f"aromatic ring N atom {atom} protonated")
    if pdb_id in ALIPHATIC_N_PROTONATION:
        atom = protonate_aliphatic_n(mol)
        changes.append(f"tertiary aliphatic N atom {atom} protonated")
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    prepared = Chem.AddHs(mol, addCoords=True)
    prepared.GetConformer().Set3D(True)
    prepared.SetProp("_Name", f"{pdb_id}_{ccd}_crystal_pH{PH}")
    prepared.SetProp("PDB_ID", pdb_id)
    prepared.SetProp("CCD_ID", ccd)
    prepared.SetProp("PROTONATION_PH", str(PH))
    prepared.SetProp("PROTONATION_NOTES", "; ".join(changes) if changes else "neutral protomer")
    charge = Chem.GetFormalCharge(prepared)
    if charge != LIGAND_CHARGES[pdb_id]:
        raise RuntimeError(f"prepared ligand charge {charge}, expected {LIGAND_CHARGES[pdb_id]}")
    new_heavy = heavy_coordinates(prepared)
    if len(original_heavy) != len(new_heavy):
        raise RuntimeError("heavy atom count changed during ligand preparation")
    rmsd = math.sqrt(sum(sum((a - b) ** 2 for a, b in zip(x, y)) for x, y in zip(original_heavy, new_heavy)) / len(new_heavy))
    if rmsd > 1e-6:
        raise RuntimeError(f"ligand heavy coordinates changed: RMSD={rmsd}")
    writer = Chem.SDWriter(str(output))
    writer.write(prepared)
    writer.close()
    return output, {
        "reused": False,
        "charge": charge,
        "atom_count": prepared.GetNumAtoms(),
        "heavy_atom_count": sum(atom.GetAtomicNum() > 1 for atom in prepared.GetAtoms()),
        "heavy_atom_rmsd_A": rmsd,
        "protonation_changes": changes or ["neutral protomer"],
    }


def validate_for_qmeku(receptor: Path, ligand: Path) -> dict[str, object]:
    from qmeku.analysis import check_residues
    from qmeku.ligand import build_complex, read_poses
    from qmeku.structure import (
        Complex, fragment_atoms, link_atoms_distance_check,
        n_electrons, pocket_atoms, residue_charge, select_residues, system_charge,
    )

    poses = read_poses(ligand)
    if len(poses) != 1:
        raise RuntimeError(f"expected one ligand pose, found {len(poses)}")
    name, molecule = poses[0]
    lig_charge = Chem.GetFormalCharge(molecule)
    complex_text = build_complex(receptor, molecule)
    with tempfile.TemporaryDirectory() as tmp_name:
        complex_path = Path(tmp_name) / "complex.pdb"
        complex_path.write_text(complex_text)
        complex_model = Complex(complex_path, "LIG")
        selected = select_residues(complex_model, 4.0)
        atoms, used = pocket_atoms(complex_model, selected)
        check_residues(complex_model, complex_model.protein_keys)
        for k in complex_model.protein_keys:
            q_res = residue_charge(complex_model, k)
            if n_electrons(fragment_atoms(complex_model, k), q_res) % 2:
                raise RuntimeError(f"odd electron count in receptor residue {complex_model.label(k)}")
        lig_atoms = complex_model.ligand_atoms()
        pocket_q = system_charge(complex_model, used)
        charge = pocket_q + lig_charge
        electrons = n_electrons(atoms + lig_atoms, charge)
        min_dist = link_atoms_distance_check(atoms + lig_atoms)
        if (not selected or n_electrons(atoms, pocket_q) % 2 or n_electrons(lig_atoms, lig_charge) % 2
                or electrons % 2 or min_dist < 0.5):
            raise RuntimeError(f"invalid QMeku model: pocket={len(selected)}, electrons={electrons}, min_dist={min_dist:.3f}")
    return {
        "pose": name,
        "ligand_charge": lig_charge,
        "pocket_residues_at_4A": len(used),
        "pocket_atoms": len(atoms),
        "model_charge": charge,
        "model_electrons": electrons,
        "min_atom_distance_A": round(min_dist, 4),
        "even_electron_model": electrons % 2 == 0,
    }


def update_manifest(rows: list[dict[str, str]], prepared_ids: set[str]) -> None:
    temporary = MANIFEST.with_suffix(".csv.part")
    for row in rows:
        if row["pdb_id"].upper() in prepared_ids:
            row["status"] = "prepared"
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(MANIFEST)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", help="optional PDB IDs; default: all downloaded benchmark entries")
    parser.add_argument("--force", action="store_true", help="overwrite existing prepared files")
    args = parser.parse_args()
    rows = read_manifest()
    wanted = {item.upper() for item in args.ids} if args.ids else None
    selected = [row for row in rows if row["pdb_id"] != "3PTB" and (wanted is None or row["pdb_id"].upper() in wanted)]
    amber_bin = ambertools_bin()
    completed: set[str] = set()
    failures: list[tuple[str, str]] = []
    for index, row in enumerate(selected, start=1):
        pdb_id = row["pdb_id"].upper()
        case = BENCHMARK / row["set"] / pdb_id
        print(f"[{index}/{len(selected)}] {pdb_id}", file=sys.stderr)
        try:
            receptor, receptor_meta = prepare_receptor(case, pdb_id, amber_bin, args.force)
            ligand, ligand_meta = prepare_ligand(case, pdb_id, row["ligand_ccd"], args.force)
            validation = validate_for_qmeku(receptor, ligand)
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
            (case / "prepared" / "preparation.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
            completed.add(pdb_id)
            print(
                f"  ready: ligand charge {validation['ligand_charge']:+d}; "
                f"pocket {validation['pocket_residues_at_4A']} residues / {validation['pocket_atoms']} atoms",
                file=sys.stderr,
            )
        except Exception as exc:
            failures.append((pdb_id, str(exc)))
            print(f"  FAILED: {exc}", file=sys.stderr)
    update_manifest(rows, completed)
    print(f"Preparation reconciliation: requested={len(selected)}, prepared={len(completed)}, failed={len(failures)}", file=sys.stderr)
    for pdb_id, message in failures:
        print(f"  {pdb_id}: {message}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
