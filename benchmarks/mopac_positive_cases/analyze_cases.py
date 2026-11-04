#!/usr/bin/env python3
"""Compile mechanistic diagnostics for positive MOPAC QMEku cases and controls."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from scipy.spatial.distance import cdist


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
ANALYSIS = ROOT / "analysis"
CONTROLS = {
    "4HBV": REPO / "benchmarks/brd4_neutral_series/pfi1_series/4HBV",
    "5F63": REPO / "benchmarks/affinity_crystal/brd4_rank/5F63",
    "3P5O": REPO / "benchmarks/affinity_crystal/core_absolute/3P5O",
    "2HB1": REPO / "benchmarks/affinity_crystal/stress/2HB1",
    "1E66": REPO / "benchmarks/affinity_crystal/stress/1E66",
}


def pdb_heavy_atoms(path: Path, records: tuple[str, ...]) -> list[tuple[np.ndarray, str, str, str]]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith(records):
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element.upper() == "H":
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        atoms.append((xyz, line[17:20].strip(), line[21:22], line[22:26].strip()))
    return atoms


def first_sdf(path: Path) -> Chem.Mol:
    mol = next((item for item in Chem.SDMolSupplier(str(path), removeHs=False) if item), None)
    if mol is None:
        raise ValueError(f"Cannot read {path}")
    return mol


def nearby_crystal_content(case: Path) -> tuple[int, str]:
    ligand = np.array([atom[0] for atom in pdb_heavy_atoms(case / "raw/ligand_crystal.pdb", ("HETATM", "ATOM  "))])
    hetero = pdb_heavy_atoms(case / "raw/structure.pdb", ("HETATM",))
    waters: set[tuple[str, str, str]] = set()
    other: dict[tuple[str, str, str], float] = {}
    for xyz, residue, chain, number in hetero:
        distance = float(cdist(ligand, [xyz]).min())
        if distance < 0.1:  # atom of the selected ligand itself
            continue
        key = (residue, chain, number)
        if residue in {"HOH", "WAT"} and distance < 4.0:
            waters.add(key)
        elif residue not in {"HOH", "WAT"} and distance < 5.0:
            other[key] = min(other.get(key, 99.0), distance)
    labels = ";".join(f"{r}:{c}{n}@{d:.2f}A" for (r, c, n), d in sorted(other.items(), key=lambda x: x[1]))
    return len(waters), labels


def diagnostics(group: str, pdb_id: str, case: Path, split_path: Path) -> dict[str, object]:
    with split_path.open(newline="") as handle:
        energy = next(csv.DictReader(handle))
    mol = first_sdf(case / "prepared/ligand.sdf")
    conformer = mol.GetConformer()
    ligand_xyz = np.array([
        [*conformer.GetAtomPosition(atom.GetIdx())]
        for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1
    ])
    receptor_xyz = np.array([atom[0] for atom in pdb_heavy_atoms(case / "prepared/receptor.pdb", ("ATOM  ", "HETATM"))])
    distances = cdist(ligand_xyz, receptor_xyz)
    waters, other_hetero = nearby_crystal_content(case)
    entry = json.loads((case / "raw/entry.json").read_text())
    title = entry.get("struct", {}).get("title", "")
    mopac_gas = float(energy["mopac_e_gas"])
    mopac_solv = float(energy["mopac_d_solv"])
    xtb_gas = float(energy["xtb_e_gas"])
    xtb_solv = float(energy["xtb_d_solv"])
    return {
        "group": group,
        "pdb_id": pdb_id,
        "target_or_structure": title,
        "ligand_charge": int(energy["ligand_charge"]),
        "pocket_charge": int(energy["pocket_charge"]),
        "pocket_residues": int(energy["pocket_residues"]),
        "pocket_atoms": int(energy["pocket_atoms"]),
        "xtb_e_gas": xtb_gas,
        "xtb_d_solv": xtb_solv,
        "xtb_e_bind": float(energy["xtb_e_bind"]),
        "mopac_e_gas": mopac_gas,
        "mopac_d_solv": mopac_solv,
        "mopac_e_bind": float(energy["mopac_e_bind"]),
        "mopac_desolv_over_abs_gas": mopac_solv / abs(mopac_gas) if mopac_gas < 0 else "",
        "mopac_minus_xtb_e_bind": float(energy["mopac_e_bind"]) - float(energy["xtb_e_bind"]),
        "mopac_minus_xtb_d_solv": mopac_solv - xtb_solv,
        "ligand_heavy_atoms": mol.GetNumHeavyAtoms(),
        "molecular_weight": Descriptors.MolWt(mol),
        "tpsa_A2": rdMolDescriptors.CalcTPSA(mol),
        "clogP": Crippen.MolLogP(mol),
        "hbond_donors": Lipinski.NumHDonors(mol),
        "hbond_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "min_heavy_distance_A": float(distances.min()),
        "heavy_contacts_below_2A": int((distances < 2.0).sum()),
        "crystal_waters_within_4A": waters,
        "other_excluded_hetero_within_5A": other_hetero,
        "source_case": case.relative_to(REPO).as_posix(),
    }


def main() -> None:
    ANALYSIS.mkdir(exist_ok=True)
    with (ROOT / "manifest.csv").open(newline="") as handle:
        positive_sources = {row["pdb_id"]: REPO / row["source_case"] for row in csv.DictReader(handle)}
    rows = []
    for pdb_id, source in sorted(positive_sources.items()):
        rows.append(diagnostics("positive", pdb_id, source, ROOT / pdb_id / "diagnostics/gas_split/poses.csv"))
    for pdb_id, source in CONTROLS.items():
        rows.append(diagnostics("matched_control", pdb_id, source, ANALYSIS / "controls" / pdb_id / "gas_split/poses.csv"))
    output = ANALYSIS / "case_diagnostics.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} diagnostic rows to {output}")


if __name__ == "__main__":
    main()
