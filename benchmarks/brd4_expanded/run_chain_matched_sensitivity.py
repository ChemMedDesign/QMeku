#!/usr/bin/env python3
"""Recalculate multichain BRD4 cases using only the ligand-matched protein chain."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent
OUTPUT = BENCHMARK / "chain_matched_sensitivity"
SHARED = BENCHMARK.parent / "affinity_crystal" / "prepare_benchmark.py"
IDS = ("4OGI", "5VBP")


def load_shared():
    spec = importlib.util.spec_from_file_location("affinity_crystal_preparation", SHARED)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {SHARED}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_chain_receptor(source: Path, destination: Path, chain: str) -> int:
    atoms: list[str] = []
    for line in source.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM  ") and line[21:22] == chain and line[16:17] in (" ", "A"):
            atoms.append(line[:16] + " " + line[17:])
    if not atoms:
        raise RuntimeError(f"no protein atoms found for chain {chain!r}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(atoms + ["TER", "END"]) + "\n", encoding="ascii")
    return len(atoms)


def one_energy(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle))


def main() -> int:
    shared = load_shared()
    amber_bin = shared.ambertools_bin()
    comparisons: list[dict[str, object]] = []
    for pdb_id in IDS:
        original = BENCHMARK / "kd_itc_series" / pdb_id
        case = OUTPUT / pdb_id
        metadata = json.loads((original / "metadata.json").read_text(encoding="utf-8"))
        chain = metadata["selected_ligand_instance"]["chain"]
        raw_receptor = case / "raw" / "receptor_protein_only_unprepared.pdb"
        source_atom_count = write_chain_receptor(original / "raw" / "structure.pdb", raw_receptor, chain)
        receptor, receptor_meta = shared.prepare_receptor(case, pdb_id, amber_bin, True)
        ligand = original / "prepared" / "ligand.sdf"
        validation = shared.validate_for_qmeku(receptor, ligand)
        outdir = case / "results"
        command = [
            str(ROOT / "bin" / "qmeku"), str(receptor), str(ligand),
            "--outdir", str(outdir), "--method", "both", "--total-only",
            "--no-gas", "--workers", "2", "--quiet",
        ]
        completed = subprocess.run(command, cwd=ROOT, check=False)
        if completed.returncode:
            raise RuntimeError(f"QMeku failed for {pdb_id} ({completed.returncode})")
        old = one_energy(original / "results" / "poses.csv")
        new = one_energy(outdir / "poses.csv")
        comparison = {
            "pdb_id": pdb_id,
            "ligand_chain": chain,
            "original_selected_chains": ";".join(metadata["selected_protein_chains_within_8A"]),
            "chain_matched_source_atoms": source_atom_count,
            "original_pocket_residues": old["pocket_residues"],
            "chain_matched_pocket_residues": new["pocket_residues"],
            "original_pocket_atoms": old["pocket_atoms"],
            "chain_matched_pocket_atoms": new["pocket_atoms"],
            "original_pocket_charge": old["pocket_charge"],
            "chain_matched_pocket_charge": new["pocket_charge"],
            "original_xtb": old["xtb_e_bind"],
            "chain_matched_xtb": new["xtb_e_bind"],
            "delta_xtb": f"{float(new['xtb_e_bind']) - float(old['xtb_e_bind']):.3f}",
            "original_mopac": old["mopac_e_bind"],
            "chain_matched_mopac": new["mopac_e_bind"],
            "delta_mopac": f"{float(new['mopac_e_bind']) - float(old['mopac_e_bind']):.3f}",
        }
        comparisons.append(comparison)
        audit = {
            "pdb_id": pdb_id,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": "chain-matched monomer sensitivity; original benchmark retained",
            "ligand_chain": chain,
            "original_selected_chains": metadata["selected_protein_chains_within_8A"],
            "receptor": receptor_meta,
            "qmeku_validation": validation,
            "qmeku_command": command,
        }
        (case / "sensitivity.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    output_csv = OUTPUT / "comparison.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    print(output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
