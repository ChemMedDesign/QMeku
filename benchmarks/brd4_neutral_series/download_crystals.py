#!/usr/bin/env python3
"""Download and organize the neutral BRD4 crystal structures into raw/ directories.

Based on the robust downloader from benchmarks/affinity_crystal/download_crystals.py.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BENCHMARK = Path(__file__).resolve().parent
MANIFEST = BENCHMARK / "manifest.csv"
USER_AGENT = "QMeku-crystal-benchmark/1.0 (RCSB targeted downloader)"

SOLVENTS = {
    "HOH", "DOD", "EDO", "DMS", "GOL", "PEG", "PGE", "PG4", "SO4", "PO4", "NO3", "CL",
    "NA", "ACT", "FMT", "MPD", "TRS", "EPE", "MES", "1PE", "P6G", "2PE", "IOD", "BR",
    "MG", "CA", "ZN", "NI", "CIT", "ACY", "IPA", "EOH", "BME", "DTT",
}


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str, destination: Path, force: bool, retries: int = 2) -> None:
    if destination.is_file() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=45) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
            if temporary.stat().st_size == 0:
                raise OSError("empty response")
            temporary.replace(destination)
            return
        except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            temporary.unlink(missing_ok=True)
            if attempt == retries:
                raise RuntimeError(f"could not download {url}: {exc}") from exc
            time.sleep(1.0 + attempt)


def first_model(lines: list[str]) -> list[str]:
    if not any(line.startswith("MODEL ") for line in lines):
        return lines
    selected: list[str] = []
    inside = False
    for line in lines:
        if line.startswith("MODEL "):
            if inside:
                break
            inside = True
            continue
        if inside and line.startswith("ENDMDL"):
            break
        if inside or line.startswith(("HEADER", "TITLE ", "COMPND", "SOURCE", "REMARK", "CONECT")):
            selected.append(line)
    selected.extend(line for line in lines if line.startswith("CONECT"))
    return selected


def choose_ligand(lines: list[str], ccd: str) -> tuple[tuple[str, str, str], list[str]]:
    groups: dict[tuple[str, str, str], list[str]] = {}
    for line in lines:
        if not line.startswith(("HETATM", "ATOM  ")):
            continue
        if line[17:20].strip() != ccd:
            continue
        key = (line[21], line[22:26].strip(), line[26].strip() or "_")
        groups.setdefault(key, []).append(line)
    if not groups:
        raise RuntimeError(f"no ligand instance found for CCD {ccd}")
    # Choose instance with most atoms or first
    chosen_key = max(groups.keys(), key=lambda k: len(groups[k]))
    return chosen_key, groups[chosen_key]


def write_clean_pdb(lines: list[str], destination: Path) -> None:
    text = "\n".join(lines).strip() + "\nEND\n"
    destination.write_text(text, encoding="ascii")


def isolate_receptor(lines: list[str], ligand_lines: list[str], destination: Path) -> list[str]:
    lig_coords = []
    for l in ligand_lines:
        if l.startswith(("ATOM", "HETATM")):
            lig_coords.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    selected = []
    chains_near = set()
    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        resname = line[17:20].strip()
        chain = line[21]
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if any((xyz[0]-lx)**2 + (xyz[1]-ly)**2 + (xyz[2]-lz)**2 < 64.0 for lx, ly, lz in lig_coords):
            chains_near.add(chain)
        selected.append(line)
    # Filter to protein chains within 8A of ligand
    final_lines = [l for l in selected if l[21] in chains_near]
    write_clean_pdb(final_lines, destination)
    return sorted(chains_near)


def build_raw_ligand_sdf(ideal_sdf: Path, crystal_pdb: Path, output_sdf: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    template = Chem.MolFromMolFile(str(ideal_sdf), removeHs=True)
    if template is None:
        raise RuntimeError(f"could not parse ideal SDF: {ideal_sdf}")
    crystal = Chem.MolFromPDBFile(str(crystal_pdb), removeHs=False, sanitize=False, proximityBonding=True)
    if crystal is None:
        raise RuntimeError(f"could not parse crystal PDB: {crystal_pdb}")
    assigned = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    Chem.SanitizeMol(assigned)
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(assigned)
    writer.close()


def process_case(row: dict[str, str], force: bool) -> None:
    pdb_id = row["pdb_id"].upper()
    ccd = row["ligand_ccd"]
    case_dir = BENCHMARK / row["set"] / pdb_id
    raw_dir = case_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    json_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    ideal_url = f"https://files.rcsb.org/ligands/download/{ccd}_ideal.sdf"

    fetch(pdb_url, raw_dir / "structure.pdb", force)
    fetch(cif_url, raw_dir / "structure.cif", force)
    fetch(json_url, raw_dir / "entry.json", force)
    fetch(ideal_url, raw_dir / f"{ccd}_ideal.sdf", force)

    pdb_lines = first_model((raw_dir / "structure.pdb").read_text(encoding="ascii", errors="replace").splitlines())
    chosen_key, lig_lines = choose_ligand(pdb_lines, ccd)
    write_clean_pdb(lig_lines, raw_dir / "ligand_crystal.pdb")

    chains = isolate_receptor(pdb_lines, lig_lines, raw_dir / "receptor_protein_only_unprepared.pdb")

    raw_sdf = raw_dir / "ligand_crystal_unprepared.sdf"
    build_raw_ligand_sdf(raw_dir / f"{ccd}_ideal.sdf", raw_dir / "ligand_crystal.pdb", raw_sdf)

    meta = {
        "pdb_id": pdb_id,
        "set": row["set"],
        "ligand_ccd": ccd,
        "selected_instance": {"chain": chosen_key[0], "resnum": chosen_key[1], "inscode": chosen_key[2]},
        "chains_within_8A": chains,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": {p.name: sha256(p) for p in raw_dir.iterdir() if p.is_file()},
    }
    (case_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"[{pdb_id}] Downloaded and raw files generated successfully (chain {chosen_key[0]}:{chosen_key[1]}).")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids", nargs="*", help="Optional PDB IDs to download")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    rows = read_manifest()
    wanted = {i.upper() for i in args.ids} if args.ids else None
    failures = 0
    for r in rows:
        if wanted and r["pdb_id"].upper() not in wanted:
            continue
        try:
            process_case(r, args.force)
        except Exception as exc:
            failures += 1
            print(f"[{r['pdb_id']}] FAILED: {exc}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
