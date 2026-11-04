#!/usr/bin/env python3
"""Download and organize the RCSB crystal structures in manifest.csv.

The script keeps deposited files under ``raw/`` and never labels a receptor as
QMeku-ready.  It extracts a ligand instance with crystal coordinates and a raw
protein-only receptor, but protonation and AMBER residue naming remain an explicit
preparation step.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


BENCHMARK = Path(__file__).resolve().parent
MANIFEST = BENCHMARK / "manifest.csv"
PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
CCD_RE = re.compile(r"^[A-Za-z0-9]{1,3}$")
USER_AGENT = "QMeku-crystal-benchmark/1.0 (RCSB targeted downloader)"


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if not PDB_RE.fullmatch(row["pdb_id"]):
            raise ValueError(f"invalid PDB ID in manifest: {row['pdb_id']!r}")
        if not CCD_RE.fullmatch(row["ligand_ccd"]):
            raise ValueError(f"invalid CCD ID in manifest: {row['ligand_ccd']!r}")
    return rows


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


def atom_coordinates(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


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
    # CONECT records are commonly after ENDMDL and are added from the full file later.
    selected.extend(line for line in lines if line.startswith("CONECT"))
    return selected


def choose_ligand(lines: list[str], ccd: str) -> tuple[tuple[str, str, str], list[str]]:
    groups: dict[tuple[str, str, str], list[str]] = {}
    for line in lines:
        if not line.startswith("HETATM") or line[17:20].strip().upper() != ccd:
            continue
        altloc = line[16:17]
        if altloc not in (" ", "A"):
            continue
        key = (line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip() or "_")
        groups.setdefault(key, []).append(line)
    if not groups:
        raise ValueError(f"CCD ligand {ccd} is absent from the selected PDB coordinates")

    def score(item: tuple[tuple[str, str, str], list[str]]) -> tuple[int, float, str]:
        key, atoms = item
        heavy = sum(1 for line in atoms if line[76:78].strip().upper() not in ("H", "D"))
        occupancy = sum(float(line[54:60].strip() or "0") for line in atoms)
        return heavy, occupancy, ":".join(key)

    key, atoms = max(groups.items(), key=score)
    return key, atoms


def write_ligand_pdb(lines: list[str], atoms: list[str], destination: Path) -> None:
    serials = {int(line[6:11]) for line in atoms}
    output = list(atoms)
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        fields = [int(value) for value in line[6:].split() if value.isdigit()]
        if fields and fields[0] in serials:
            kept = [value for value in fields if value in serials]
            if len(kept) > 1:
                output.append("CONECT" + "".join(f"{value:5d}" for value in kept))
    output.append("END")
    destination.write_text("\n".join(output) + "\n", encoding="ascii")


def contacting_chains(lines: list[str], ligand_atoms: list[str], cutoff: float = 8.0) -> list[str]:
    ligand_xyz = [
        atom_coordinates(line)
        for line in ligand_atoms
        if line[76:78].strip().upper() not in ("H", "D")
    ]
    cutoff2 = cutoff * cutoff
    chains: set[str] = set()
    for line in lines:
        if not line.startswith("ATOM  ") or line[16:17] not in (" ", "A"):
            continue
        xyz = atom_coordinates(line)
        if any(sum((a - b) ** 2 for a, b in zip(xyz, lig)) <= cutoff2 for lig in ligand_xyz):
            chains.add(line[21:22])
    if not chains:
        chains = {line[21:22] for line in lines if line.startswith("ATOM  ")}
    return sorted(chains)


def write_raw_receptor(lines: list[str], chains: list[str], destination: Path) -> None:
    by_chain: dict[str, list[str]] = {chain: [] for chain in chains}
    for line in lines:
        if line.startswith("ATOM  ") and line[21:22] in by_chain and line[16:17] in (" ", "A"):
            by_chain[line[21:22]].append(line)
    output: list[str] = []
    for chain in chains:
        output.extend(by_chain[chain])
        if by_chain[chain]:
            output.append("TER")
    output.append("END")
    destination.write_text("\n".join(output) + "\n", encoding="ascii")


def convert_crystal_ligand(source: Path, ideal_sdf: Path, destination: Path) -> str:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        template = Chem.MolFromMolFile(str(ideal_sdf), removeHs=True)
        crystal = Chem.MolFromPDBFile(str(source), removeHs=False, sanitize=False, proximityBonding=True)
        if template is None or crystal is None:
            raise ValueError("RDKit could not read the CCD template or crystal ligand")
        assigned = AllChem.AssignBondOrdersFromTemplate(template, crystal)
        Chem.SanitizeMol(assigned)
        writer = Chem.SDWriter(str(destination))
        writer.write(assigned)
        writer.close()
        if destination.is_file() and destination.stat().st_size > 0:
            return "created_with_RDKIT_CCD_template_and_crystal_coordinates; verify protonation"
    except Exception as exc:
        rdkit_error = str(exc).replace("\n", " ")[:200]
        destination.unlink(missing_ok=True)
    else:
        rdkit_error = "empty RDKit output"

    executable = shutil.which("obabel")
    if executable is None:
        return f"not_created: RDKit failed ({rdkit_error}); obabel is not installed"
    completed = subprocess.run(
        [executable, "-ipdb", str(source), "-osdf", "-O", str(destination)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode or not destination.is_file() or destination.stat().st_size == 0:
        destination.unlink(missing_ok=True)
        message = completed.stderr.strip().replace("\n", " ")[:300]
        return f"not_created: RDKit failed ({rdkit_error}); obabel failed ({message})"
    return (
        "created_with_obabel_from_crystal_coordinates_after_RDKIT_failure; "
        f"verify bond orders, formal charge, and protonation (RDKit: {rdkit_error})"
    )


def urls_for(row: dict[str, str]) -> dict[str, str]:
    pdb_id = row["pdb_id"].upper()
    ccd = row["ligand_ccd"].upper()
    lower = pdb_id.lower()
    urls = {
        "entry_json": f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
        "pdb": f"https://files.rcsb.org/download/{pdb_id}.pdb",
        "cif": f"https://files.rcsb.org/download/{pdb_id}.cif",
        "ligand_ideal_sdf": f"https://files.rcsb.org/ligands/download/{ccd}_ideal.sdf",
        "validation_pdf": (
            f"https://files.rcsb.org/pub/pdb/validation_reports/{lower[1:3]}/"
            f"{lower}/{lower}_full_validation.pdf"
        ),
    }
    assembly_id = row.get("assembly_id", "").strip()
    if assembly_id:
        urls["assembly_pdb_gz"] = f"https://files.rcsb.org/download/{pdb_id}.pdb{assembly_id}.gz"
    return urls


def download_one(row: dict[str, str], force: bool, validation: bool) -> dict[str, object]:
    pdb_id = row["pdb_id"].upper()
    ccd = row["ligand_ccd"].upper()
    case = BENCHMARK / row["set"] / pdb_id
    raw = case / "raw"
    prepared = case / "prepared"
    reports = case / "validation"
    for directory in (raw, prepared, reports, case / "results"):
        directory.mkdir(parents=True, exist_ok=True)

    urls = urls_for(row)
    targets = {
        "entry_json": raw / "entry.json",
        "pdb": raw / "structure.pdb",
        "cif": raw / "structure.cif",
        "ligand_ideal_sdf": raw / f"{ccd}_ideal.sdf",
    }
    if validation:
        targets["validation_pdf"] = reports / "wwpdb_validation.pdf"
    validation_warning = None
    for name, destination in targets.items():
        try:
            fetch(urls[name], destination, force)
        except RuntimeError as exc:
            if name != "validation_pdf":
                raise
            validation_warning = str(exc)

    coordinate_pdb = targets["pdb"]
    assembly_id = row.get("assembly_id", "").strip()
    if assembly_id:
        compressed = raw / f"assembly{assembly_id}.pdb.gz"
        expanded = raw / f"assembly{assembly_id}.pdb"
        fetch(urls["assembly_pdb_gz"], compressed, force)
        if force or not expanded.is_file():
            with gzip.open(compressed, "rb") as source, expanded.open("wb") as output:
                shutil.copyfileobj(source, output)
        coordinate_pdb = expanded

    all_lines = coordinate_pdb.read_text(encoding="ascii", errors="replace").splitlines()
    lines = first_model(all_lines)
    instance, ligand_atoms = choose_ligand(lines, ccd)
    ligand_pdb = raw / "ligand_crystal.pdb"
    write_ligand_pdb(lines, ligand_atoms, ligand_pdb)
    chains = contacting_chains(lines, ligand_atoms)
    receptor_pdb = raw / "receptor_protein_only_unprepared.pdb"
    write_raw_receptor(lines, chains, receptor_pdb)
    ligand_sdf = raw / "ligand_crystal_unprepared.sdf"
    ligand_conversion = convert_crystal_ligand(ligand_pdb, targets["ligand_ideal_sdf"], ligand_sdf)

    files = [path for path in (*targets.values(), coordinate_pdb, ligand_pdb, receptor_pdb, ligand_sdf) if path.is_file()]
    metadata = {
        "pdb_id": pdb_id,
        "ligand_ccd": ccd,
        "set": row["set"],
        "accessed_utc": datetime.now(timezone.utc).isoformat(),
        "source": "RCSB Protein Data Bank",
        "urls": urls,
        "coordinate_source": coordinate_pdb.relative_to(case).as_posix(),
        "assembly_id": assembly_id or None,
        "selected_ligand_instance": {"chain": instance[0], "residue_number": instance[1], "insertion_code": instance[2]},
        "selected_protein_chains_within_8A": [chain or "_" for chain in chains],
        "ligand_atom_count": len(ligand_atoms),
        "ligand_sdf_status": ligand_conversion,
        "validation_warning": validation_warning,
        "preparation_status": "raw_only; protonation, AMBER residue naming, and ligand chemistry require validation",
        "sha256": {path.relative_to(case).as_posix(): sha256(path) for path in sorted(set(files))},
    }
    (case / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def update_manifest_status(rows: list[dict[str, str]], completed: set[str]) -> None:
    if not completed:
        return
    fieldnames = list(rows[0])
    temporary = MANIFEST.with_suffix(".csv.part")
    for row in rows:
        if row["pdb_id"].upper() in completed and row["status"] not in ("prepared", "completed"):
            row["status"] = "downloaded_raw"
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(MANIFEST)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ids", nargs="*", help="optional PDB IDs; default: every non-completed manifest row")
    parser.add_argument("--force", action="store_true", help="redownload and overwrite raw database files")
    parser.add_argument("--no-validation", action="store_true", help="skip wwPDB validation PDF downloads")
    parser.add_argument("--dry-run", action="store_true", help="print destinations and URLs without downloading")
    args = parser.parse_args()

    rows = read_manifest()
    wanted = {item.upper() for item in args.ids} if args.ids else None
    selected = [
        row for row in rows
        if (wanted is None or row["pdb_id"].upper() in wanted)
        and (args.force or row["status"] not in ("downloaded_raw", "prepared", "completed"))
    ]
    if wanted is not None:
        known = {row["pdb_id"].upper() for row in rows}
        unknown = sorted(wanted - known)
        if unknown:
            parser.error("IDs absent from manifest: " + ", ".join(unknown))

    if args.dry_run:
        for row in selected:
            case = BENCHMARK / row["set"] / row["pdb_id"].upper()
            print(f"{row['pdb_id'].upper()} -> {case}")
            for name, url in urls_for(row).items():
                if name != "validation_pdf" or not args.no_validation:
                    print(f"  {name}: {url}")
        print(f"Planned targeted retrievals: {len(selected)} PDB entries")
        return 0

    completed: set[str] = set()
    failures: list[tuple[str, str]] = []
    for index, row in enumerate(selected, start=1):
        pdb_id = row["pdb_id"].upper()
        print(f"[{index}/{len(selected)}] {pdb_id}", file=sys.stderr)
        try:
            metadata = download_one(row, args.force, not args.no_validation)
            completed.add(pdb_id)
            instance = metadata["selected_ligand_instance"]
            print(
                f"  ligand {row['ligand_ccd']} {instance['chain']}:{instance['residue_number']}; "
                f"raw download complete",
                file=sys.stderr,
            )
        except Exception as exc:  # keep the bounded batch moving and report every failed ID
            failures.append((pdb_id, str(exc)))
            print(f"  FAILED: {exc}", file=sys.stderr)
    update_manifest_status(rows, completed)

    print(
        f"Count reconciliation: requested={len(selected)}, downloaded={len(completed)}, failed={len(failures)}",
        file=sys.stderr,
    )
    for pdb_id, message in failures:
        print(f"  {pdb_id}: {message}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
