#!/usr/bin/env python3
"""Run QMeku consistently over prepared crystal benchmark complexes."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = Path(__file__).resolve().parent


def read_manifest() -> list[dict[str, str]]:
    with (BENCHMARK / "manifest.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*", help="optional PDB IDs to run")
    parser.add_argument("--method", choices=("xtb", "mopac", "both"), default="both")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    wanted = {value.upper() for value in args.ids} if args.ids else None
    failures = 0
    skipped = 0
    for row in read_manifest():
        pdb_id = row["pdb_id"].upper()
        if wanted is not None and pdb_id not in wanted:
            continue
        receptor = ROOT / row["receptor_path"]
        ligand = ROOT / row["ligand_path"]
        outdir = BENCHMARK / row["set"] / pdb_id / "results"
        result_csv = outdir / "poses.csv"
        stale = result_csv.is_file() and any(
            p.is_file() and p.stat().st_mtime > result_csv.stat().st_mtime
            for p in (receptor, ligand)
        )
        if result_csv.exists() and not stale and not args.force:
            print(f"{pdb_id}: already completed", file=sys.stderr)
            continue
        if not receptor.is_file() or not ligand.is_file():
            print(f"{pdb_id}: skipped; prepared receptor or ligand is missing", file=sys.stderr)
            skipped += 1
            continue
        command = [
            str(ROOT / "bin" / "qmeku"),
            str(receptor),
            str(ligand),
            "--outdir", str(outdir),
            "--method", args.method,
            "--total-only",
            "--no-gas",
            "--workers", str(args.workers),
            "--quiet",
        ]
        print(f"{pdb_id}: running", file=sys.stderr)
        completed = subprocess.run(command, cwd=ROOT, check=False)
        if completed.returncode:
            failures += 1
            print(f"{pdb_id}: failed with exit code {completed.returncode}", file=sys.stderr)

    print(f"Finished: failures={failures}, missing_inputs={skipped}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
