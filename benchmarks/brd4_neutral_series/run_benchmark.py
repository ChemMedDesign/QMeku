#!/usr/bin/env python3
"""Run the current QMeku over prepared neutral BRD4 crystals."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", help="optional PDB IDs")
    parser.add_argument("--method", choices=("xtb", "mopac", "both"), default="both")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    wanted = {value.upper() for value in args.ids} if args.ids else None
    with (BENCHMARK / "manifest.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    failures = 0
    missing = 0
    for row in rows:
        pdb_id = row["pdb_id"].upper()
        if wanted is not None and pdb_id not in wanted:
            continue
        receptor = ROOT / row["receptor_path"]
        ligand = ROOT / row["ligand_path"]
        outdir = BENCHMARK / row["set"] / pdb_id / "results"
        result = outdir / "poses.csv"
        stale = result.is_file() and any(
            path.is_file() and path.stat().st_mtime > result.stat().st_mtime
            for path in (receptor, ligand)
        )
        if result.is_file() and not stale and not args.force:
            print(f"{pdb_id}: already completed", file=sys.stderr)
            continue
        if not receptor.is_file() or not ligand.is_file():
            print(f"{pdb_id}: missing prepared input", file=sys.stderr)
            missing += 1
            continue
        command = [
            str(ROOT / "bin" / "qmeku"), str(receptor), str(ligand),
            "--outdir", str(outdir), "--method", args.method,
            "--total-only", "--no-gas", "--workers", str(args.workers), "--quiet",
        ]
        print(f"{pdb_id}: running {args.method}", file=sys.stderr)
        completed = subprocess.run(command, cwd=ROOT, check=False)
        if completed.returncode:
            failures += 1
            print(f"{pdb_id}: failed ({completed.returncode})", file=sys.stderr)
    print(f"Run reconciliation: failures={failures}, missing_inputs={missing}", file=sys.stderr)
    return 1 if failures or missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
