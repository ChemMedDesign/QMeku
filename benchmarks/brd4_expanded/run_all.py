#!/usr/bin/env python3
"""Resume the complete BRD4 benchmark pipeline, stopping visibly on any stage."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


BENCHMARK = Path(__file__).resolve().parent
ROOT = BENCHMARK.parents[1]


def run_stage(name: str, command: list[str], status: dict[str, object]) -> None:
    started = datetime.now(timezone.utc).isoformat()
    completed = subprocess.run(command, cwd=ROOT, check=False)
    status["stages"].append({
        "name": name,
        "command": command,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "returncode": completed.returncode,
    })
    (BENCHMARK / "pipeline_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if completed.returncode:
        raise SystemExit(f"pipeline stopped at {name} (exit {completed.returncode})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=("xtb", "mopac", "both"), default="both")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--force-prepare", action="store_true")
    parser.add_argument("--force-run", action="store_true")
    args = parser.parse_args()
    python = sys.executable
    status: dict[str, object] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "stages": [],
    }
    download = [python, str(BENCHMARK / "download_crystals.py"), "--no-validation"]
    prepare = [python, str(BENCHMARK / "prepare_benchmark.py")]
    calculate = [
        python, str(BENCHMARK / "run_benchmark.py"),
        "--method", args.method, "--workers", str(args.workers),
    ]
    if args.force_download:
        download.append("--force")
    if args.force_prepare:
        prepare.append("--force")
    if args.force_run:
        calculate.append("--force")
    run_stage("download", download, status)
    run_stage("prepare", prepare, status)
    run_stage("calculate", calculate, status)
    run_stage("collect", [python, str(BENCHMARK / "collect_results.py")], status)
    status["completed_utc"] = datetime.now(timezone.utc).isoformat()
    (BENCHMARK / "pipeline_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
