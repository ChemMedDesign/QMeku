#!/usr/bin/env python3
"""Download the bounded BRD4 benchmark manifest from official RCSB endpoints."""

from __future__ import annotations

import importlib.util
from pathlib import Path


BENCHMARK = Path(__file__).resolve().parent
SHARED = BENCHMARK.parent / "affinity_crystal" / "download_crystals.py"


def main() -> int:
    spec = importlib.util.spec_from_file_location("affinity_crystal_downloader", SHARED)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load shared downloader: {SHARED}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.BENCHMARK = BENCHMARK
    module.MANIFEST = BENCHMARK / "manifest.csv"
    module.USER_AGENT = "QMeku-BRD4-expanded-benchmark/1.0 (RCSB targeted downloader)"
    return module.main()


if __name__ == "__main__":
    raise SystemExit(main())
