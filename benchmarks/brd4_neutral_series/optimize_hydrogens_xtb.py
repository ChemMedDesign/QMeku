#!/usr/bin/env python3
"""Optimize only hydrogens of a QMEku binding-site cluster and calculate its Hessian.

The receptor pocket and crystallographic ligand are combined into one cluster.
All heavy atoms are exactly fixed with xTB's $fix directive; hydrogens are free.
This estimates local vibrational modes only, not a complete binding free energy.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qmeku import structure as pk  # noqa: E402


BENCHMARK = Path(__file__).resolve().parent
MANIFEST = BENCHMARK / "manifest.csv"


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_ligand(path: Path) -> tuple[list[tuple[str, tuple[float, float, float], str]], int]:
    mol = next((m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None), None)
    if mol is None or mol.GetNumConformers() != 1:
        raise ValueError(f"Could not read one ligand conformer from {path}")
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append((atom.GetSymbol(), (pos.x, pos.y, pos.z), f"ligand:{atom.GetIdx()+1}"))
    return atoms, Chem.GetFormalCharge(mol)


def build_qmeku_cluster(receptor_path: Path, ligand_path: Path, cutoff: float, workdir: Path):
    ligand, ligand_charge = load_ligand(ligand_path)
    combined = workdir / "receptor_ligand_for_qmeku.pdb"
    lines = [line for line in receptor_path.read_text().splitlines() if line.startswith(("ATOM  ", "HETATM"))]
    serial = len(lines)
    for index, (element, xyz, _) in enumerate(ligand, 1):
        serial += 1
        atom_name = f"{element}{index}"[:4]
        lines.append(
            f"HETATM{serial:5d} {atom_name:>4s} LIG Z{1:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element.upper():>2s}"
        )
    combined.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    receptor = pk.Complex(str(combined), ligand_resname="LIG")
    selected = pk.select_residues(receptor, cutoff=cutoff)
    pocket, used = pk.pocket_atoms(receptor, selected)
    charge = pk.system_charge(receptor, used)
    atoms = pocket + ligand
    total_charge = charge + ligand_charge
    if pk.n_electrons(atoms, total_charge) % 2:
        raise ValueError("Cluster has an odd number of electrons; verify protonation and charge")
    return atoms, total_charge, len(used), charge, ligand_charge


def write_xyz(path: Path, atoms: list[tuple[str, tuple[float, float, float], str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(atoms)}\nQMEku binding-site cluster; heavy atoms fixed, H optimized\n")
        for element, xyz, _ in atoms:
            handle.write(f"{element:2s} {xyz[0]: .10f} {xyz[1]: .10f} {xyz[2]: .10f}\n")


def write_xcontrol(path: Path, atoms: list[tuple[str, tuple[float, float, float], str]]) -> list[int]:
    fixed = [i for i, (element, _, _) in enumerate(atoms, 1) if element.upper() != "H"]
    ranges: list[str] = []
    start = previous = fixed[0]
    for value in fixed[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    path.write_text("$fix\n  atoms: " + ",".join(ranges) + "\n$end\n", encoding="utf-8")
    return fixed


def run_xtb(executable: str, xyz: Path, outdir: Path, charge: int, control: Path, task: str) -> str:
    cmd = [executable, xyz.name, "--gfn", "2", "--chrg", str(charge), "--input", control.name]
    if task == "opt":
        cmd += ["--opt", "normal"]
    elif task == "hess":
        cmd += ["--hess"]
    else:
        raise ValueError(task)
    proc = subprocess.run(cmd, cwd=outdir, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (outdir / f"xtb_{task}.log").write_text(proc.stdout, encoding="utf-8")
    if proc.returncode != 0 or "normal termination of xtb" not in proc.stdout.lower():
        raise RuntimeError(f"xTB {task} failed; inspect {outdir / f'xtb_{task}.log'}")
    return proc.stdout


def read_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2+n]:
        fields = line.split()
        atoms.append((fields[0], np.array([float(x) for x in fields[1:4]])))
    if len(atoms) != n:
        raise ValueError(f"Malformed XYZ: {path}")
    return atoms


def find_optimized_xyz(outdir: Path) -> Path:
    candidates = [outdir / "xtbopt.xyz", outdir / "xtbopt.log"]
    for candidate in candidates:
        if candidate.name.endswith(".xyz") and candidate.is_file():
            return candidate
    raise FileNotFoundError("xTB optimization did not produce xtbopt.xyz")


def run_case(row: dict[str, str], args) -> dict[str, object]:
    pdb_id = row["pdb_id"].upper()
    case = BENCHMARK / row["set"] / pdb_id
    receptor_path = case / "prepared/receptor.pdb"
    ligand_path = case / "prepared/ligand.sdf"
    outdir = Path(args.outdir).resolve() / pdb_id
    outdir.mkdir(parents=True, exist_ok=True)
    atoms, charge, n_residues, pocket_charge, ligand_charge = build_qmeku_cluster(
        receptor_path, ligand_path, args.cutoff, outdir
    )
    input_xyz = outdir / "cluster_initial.xyz"
    control = outdir / "fix_heavy.xcontrol"
    write_xyz(input_xyz, atoms)
    fixed = write_xcontrol(control, atoms)
    (outdir / "atom_map.csv").write_text(
        "index,element,label,fixed\n" + "".join(
            f'{i},{element},{label},{str(element.upper() != "H").lower()}\n'
            for i, (element, _, label) in enumerate(atoms, 1)
        ), encoding="utf-8"
    )
    (outdir / "cluster.chrg").write_text(f"{charge}\n", encoding="utf-8")

    exe = shutil.which(args.xtb)
    if exe is None:
        raise FileNotFoundError(f"xTB executable not found: {args.xtb}")
    opt_output = run_xtb(exe, input_xyz, outdir, charge, control, "opt")
    optimized = find_optimized_xyz(outdir)
    initial_coords = np.array([xyz for _, xyz, _ in atoms])
    opt_atoms = read_xyz(optimized)
    if len(opt_atoms) != len(atoms):
        raise ValueError("Optimized geometry atom count differs from input")
    opt_coords = np.array([xyz for _, xyz in opt_atoms])
    displacement = np.linalg.norm(opt_coords - initial_coords, axis=1)
    max_fixed = float(displacement[np.array(fixed)-1].max(initial=0.0))
    max_hydrogen = float(displacement[[i for i,(el,_,_) in enumerate(atoms) if el.upper()=="H"]].max(initial=0.0))
    if max_fixed > 1e-5:
        raise RuntimeError(f"xTB moved fixed heavy atoms by {max_fixed:.3g} Å")

    # Use optimized coordinates as Hessian input. Keep the same exact fixing so the Hessian
    # contains only the movable hydrogen subspace in xTB versions that honor $fix for hessians.
    hess_xyz = outdir / "cluster_optimized.xyz"
    hess_xyz.write_text(optimized.read_text(encoding="utf-8"), encoding="utf-8")
    hess_output = run_xtb(exe, hess_xyz, outdir, charge, control, "hess")
    hessian_path = outdir / "hessian"
    if not hessian_path.is_file():
        hessian_path = outdir / "hessian.out"
    summary = {
        "pdb_id": pdb_id,
        "method": "GFN2-xTB",
        "temperature_K": args.temperature,
        "cutoff_A": args.cutoff,
        "cluster_atoms": len(atoms),
        "cluster_heavy_atoms_fixed": len(fixed),
        "cluster_hydrogens_optimized": len(atoms)-len(fixed),
        "pocket_residues": n_residues,
        "pocket_charge": pocket_charge,
        "ligand_charge": ligand_charge,
        "cluster_charge": charge,
        "max_fixed_heavy_displacement_A": max_fixed,
        "max_hydrogen_displacement_A": max_hydrogen,
        "optimized_xyz": str(hess_xyz),
        "hessian_file": str(hessian_path) if hessian_path.is_file() else None,
        "optimization_log": str(outdir / "xtb_opt.log"),
        "hessian_log": str(outdir / "xtb_hess.log"),
        "optimization_normal_termination": "normal termination of xtb" in opt_output.lower(),
        "hessian_normal_termination": "normal termination of xtb" in hess_output.lower(),
        "thermochemistry_note": "This script computes the constrained cluster Hessian only; derive thermal corrections from its vibrational modes. It is not a binding free energy.",
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="+", default=["4E96", "4NUD", "6FO5"])
    parser.add_argument("--cutoff", type=float, default=4.0)
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--xtb", default="xtb")
    parser.add_argument("--outdir", default=str(BENCHMARK / "hydrogen_only_xtb"))
    args = parser.parse_args()
    wanted = {pdb.upper() for pdb in args.ids}
    rows = {row["pdb_id"].upper(): row for row in read_manifest()}
    missing = wanted - rows.keys()
    if missing:
        parser.error("PDB IDs absent from manifest: " + ", ".join(sorted(missing)))
    summaries = []
    for pdb_id in sorted(wanted):
        print(f"[{pdb_id}] build pocket and optimize H atoms", file=sys.stderr, flush=True)
        try:
            summary = run_case(rows[pdb_id], args)
            summaries.append(summary)
            print(
                f"  atoms={summary['cluster_atoms']} fixed-heavy={summary['cluster_heavy_atoms_fixed']} "
                f"mobile-H={summary['cluster_hydrogens_optimized']} "
                f"max-heavy-shift={summary['max_fixed_heavy_displacement_A']:.2e} Å",
                file=sys.stderr, flush=True,
            )
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr, flush=True)
            return 1
    output = Path(args.outdir).resolve() / "summary.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
