# QMeku: approximate binding energy and per-residue decomposition of docking poses.
# Copyright (C) 2026 SB2CS lab
# SPDX-License-Identifier: GPL-3.0-or-later
# Free software under the GNU General Public License, either version 3 or (at
# your option) any later version; see the LICENSE file.
# Distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
"""Command line of QMeku:  qmeku RECEPTOR.pdb POSES.sdf --outdir RESULTS"""
import argparse
import shlex
import sys
import tempfile
import time
from pathlib import Path

from . import __version__, engines
from .analysis import analyze_ligand, pocket_residues
from .ligand import QmekuError, prepare_receptor, read_poses
from .report import format_log, summarize_residues, write_outputs


BANNER = r"""
  +----------------------------------------------------------------------------+
  |                                                                            |
  |    ____  __  __        _                                                   |
  |   / __ \|  \/  | ___  | | ___   _                                          |
  |  | |  | | |\/| |/ _ \ | |/ / | | |        QMeku v.{version:<6}                   |
  |  | |__| | |  | |  __/ |   <| |_| |        SB2CS lab                        |
  |   \___\_\_|  |_|\___| |_|\_\\__,_|                                         |
  |                                                                            |
  |   Semiempirical binding energy decomposition of docking poses              |
  |   GFN2-xTB (ALPB) and PM6-D3H4X (COSMO)                                    |
  +----------------------------------------------------------------------------+"""

EXAMPLES = """\
QMeku examples
=================

1. One pose (xtb + MOPAC, water, with the gas / solvation balance):
     qmeku examples/3ptb/receptor.pdb examples/3ptb/crystal.sdf --outdir out

2. Many poses, xtb only, no gas phase (faster), per-residue average over the poses with E_bind <= 0:
     qmeku receptor.pdb poses.sdf --outdir out --method xtb --no-gas --summary-max-ebind 0

3. Only MOPAC, gas phase (no implicit solvent):
     qmeku receptor.pdb poses.sdf --outdir out --method mopac --solvent none

4. Wider pocket (5 A), each pose with its own residues, only E_bind (no per-residue decomposition):
     qmeku receptor.pdb poses.sdf --outdir out --cutoff 5.0 --pocket per-pose --total-only

5. Protonate the receptor with PDB2PQR first, charged ligand, 8 simultaneous calculations:
     qmeku receptor_raw.pdb poses.sdf --outdir out --prepare --ligand-charge -1 --workers 8

6. Write the log in Spanish (English is the default):
     qmeku receptor.pdb poses.sdf --outdir out --lang es

Inputs:  RECEPTOR.pdb (hydrogens, AMBER residue names) and POSES.sdf (3D poses, one molecule = one pose).
Outputs: CSV files and qmeku.log in --outdir. Run 'qmeku --help' for every option.
"""


def command_line(argv):
    """The command as typed, quoted so it can be pasted back into a shell."""
    return " ".join(["qmeku"] + [shlex.quote(a) for a in argv])


def build_parser():
    ap = argparse.ArgumentParser(
        prog="qmeku",
        description="QMeku: approximate binding energy (E_bind) of docking poses with xtb and/or MOPAC, with the "
                    "per-residue decomposition, the closest distance to the ligand, the gas / solvation balance and "
                    "the GFN2 energy components (electrostatics, dispersion, repulsion, XC, solvation).",
        epilog="Energies are single points without geometry optimisation, in kcal/mol; they are not free energies.")
    ap.add_argument("--examples", action="store_true", help="show usage examples and exit")
    ap.add_argument("receptor", help="receptor PDB with hydrogens and AMBER residue names (HID/HIE/HIP, CYX); "
                                     "only amino acids (use --prepare to protonate with PDB2PQR)")
    ap.add_argument("poses", help="SDF/MOL file with the ligand poses in 3D (one molecule = one pose; "
                                  "every conformer of a molecule is a pose)")
    ap.add_argument("--outdir", required=True, help="directory for the CSV files and qmeku.log")
    ap.add_argument("--method", choices=("xtb", "mopac", "both"), default="both",
                    help="semiempirical program: xtb = GFN2-xTB (ALPB), mopac = PM6-D3H4X (COSMO) (default both)")
    ap.add_argument("--solvent", default="water", help="implicit solvent: 'water' (default) or 'none' (gas phase)")
    ap.add_argument("--cutoff", type=float, default=4.0,
                    help="a residue is in the pocket if any atom is within this distance (A) of the ligand (default 4.0)")
    ap.add_argument("--pocket", choices=("common", "per-pose"), default=None,
                    help="common = the same residues for every pose (the union of those within the cutoff in any pose; "
                         "default with several poses; needed to compare energies between poses); per-pose = each pose its own")
    ap.add_argument("--no-gas", action="store_true",
                    help="skip the gas-phase calculations (faster; no gas / solvation balance)")
    ap.add_argument("--ligand-charge", type=int, default=None, help="net charge of the ligand (default: RDKit formal charge)")
    ap.add_argument("--workers", type=int, default=None, help="simultaneous single points (default: all CPUs)")
    ap.add_argument("--no-concap", action="store_true",
                    help="do not subtract the peptide groups shared by neighbouring residue fragments (MFCC concaps); "
                         "by default they are subtracted once so they are not counted in two residues")
    ap.add_argument("--total-only", action="store_true", help="only E_bind, no per-residue decomposition")
    ap.add_argument("--summary-max-ebind", type=float, default=None,
                    help="leave out of the per-residue average the poses whose E_bind (first method) is above this value, "
                         "e.g. 0: clashing poses distort averages")
    ap.add_argument("--no-pymol", action="store_true",
                    help="do not write the PyMOL view of the per-residue decomposition (pymol/ folder: .pml scripts, and "
                         ".png when PyMOL is installed). By default it is written whenever the decomposition is done")
    ap.add_argument("--pymol-poses", type=int, default=3,
                    help="number of best poses (lowest E_bind) drawn in PyMOL (default 3; 0 = all)")
    ap.add_argument("--pymol-method", choices=("xtb", "mopac", "both"), default="both",
                    help="method(s) drawn in PyMOL, among those computed (default both: one scene per method)")
    ap.add_argument("--top", type=int, default=20, help="residues shown in the log tables (default 20; 0 = all)")
    ap.add_argument("--lang", choices=("en", "es"), default="en", help="language of the log (default en)")
    ap.add_argument("--prepare", action="store_true", help="protonate the receptor with PDB2PQR first (needs pdb2pqr)")
    ap.add_argument("--quiet", action="store_true", help="do not print the log (it is always written to qmeku.log)")
    ap.add_argument("--version", action="version", version=f"qmeku {__version__}")
    ap.add_argument("--timeout", type=float, default=600,
                    help="timeout in seconds per xtb/mopac/pdb2pqr call (default 600; 0 = no limit)")
    return ap


def run(args):
    if args.timeout is not None and args.timeout < 0:
        raise QmekuError(f"--timeout must be >= 0 (got {args.timeout})")
    timeout = None if args.timeout == 0 else args.timeout
    avail = engines.available()
    wanted = ["xtb", "mopac"] if args.method == "both" else [args.method]
    methods = [m for m in wanted if avail[m]]
    for m in wanted:
        if not avail[m]:
            msg = f"{m} is not installed (set {m.upper()}_BIN or put it in the PATH)"
            if args.method == "both" and methods:
                print("warning: " + msg + "; continuing without it", file=sys.stderr)
            else:
                raise QmekuError(msg)
    solvent_on = args.solvent not in ("none", "gas")
    if solvent_on and "mopac" in methods and args.solvent != "water":
        raise QmekuError(f"MOPAC solvent must be 'water' (got {args.solvent!r})")
    versions = engines.versions()

    with tempfile.TemporaryDirectory() as tmp:
        receptor = args.receptor
        if args.prepare:
            receptor = prepare_receptor(receptor, Path(tmp) / "receptor_h.pdb", timeout=timeout)
        receptor_text = Path(receptor).read_text()          # kept for the PyMOL files (the prepared one is temporary)
        poses = read_poses(args.poses)
        pose_mols = {}
        n_poses = sum(m.GetNumConformers() for _, m in poses)
        mode = args.pocket or ("common" if n_poses > 1 else "per-pose")
        common = None
        if mode == "common":
            common = []
            for _, m in poses:
                common += [k for k in pocket_residues(receptor, m, args.cutoff) if k not in common]
        results, failed, done = [], [], 0
        for name, mol in poses:
            try:
                for pose_name, res in analyze_ligand(
                        receptor, mol, ligand_charge=args.ligand_charge, name=name,
                        methods=tuple(methods), solvent=args.solvent, gas=not args.no_gas, cutoff=args.cutoff,
                        pairs=not args.total_only, workers=args.workers, residue_keys=common,
                        concap=not args.no_concap, timeout=timeout):
                    results.append((pose_name, res))
                    cid = (mol.GetConformer().GetId() if pose_name == name
                           else int(pose_name.rsplit("_c", 1)[1]))
                    pose_mols[pose_name] = (mol, cid)
                    done += 1
                    e = ", ".join(f"{m} {res['methods'][m]['e_bind']:.2f}" for m in methods)
                    print(f"[{done}/{n_poses}] {pose_name}: E_bind {e} kcal/mol ({res['seconds']:.0f} s)", file=sys.stderr)
            except Exception as exc:       # report and go on with the next pose
                failed.append((name, str(exc)[:300]))
                print(f"{name}: failed: {str(exc)[:200]}", file=sys.stderr)
    if not results:
        detail = f" ({failed[0][0]}: {failed[0][1]})" if failed else ""
        raise QmekuError("all poses failed" + detail)
    meta = {"version": __version__, "receptor": str(args.receptor), "poses": str(args.poses), "versions": versions,
            "solvent": args.solvent, "solvent_on": solvent_on,
            "pocket": mode, "cutoff": args.cutoff, "command": args.command_line}
    log = format_log(results, methods, meta, lang=args.lang, top=args.top, max_ebind=args.summary_max_ebind, failed=failed)
    header = BANNER.format(version=__version__) + f"\n  Command: {args.command_line}\n\n"
    write_outputs(args.outdir, results, methods, args.summary_max_ebind, log_text=header + log)
    if not args.quiet:
        print(log)
    if not args.total_only and not args.no_pymol:
        try:                                    # the PyMOL view must never break the calculation
            from .pymol_view import write_scenes
            pm = [m for m in methods if args.pymol_method in ("both", m)]
            if not pm:
                raise ValueError(f"--pymol-method {args.pymol_method} was not computed (methods run: {', '.join(methods)})")
            summ = summarize_residues(results, methods, args.summary_max_ebind) if len(results) > 1 else None
            for msg in write_scenes(args.outdir, receptor_text, results, pm, pose_mols, summary=summ,
                                    n_poses=args.pymol_poses, lang=args.lang):
                print(msg, file=sys.stderr)
        except Exception as exc:
            print(f"warning: PyMOL view not written ({str(exc)[:150]})", file=sys.stderr)
    print(f"{'Archivos escritos en' if args.lang == 'es' else 'Files written to'} {args.outdir}/", file=sys.stderr)
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if "--examples" in argv:
        print(BANNER.format(version=__version__))
        print(EXAMPLES)
        return 0
    args = build_parser().parse_args(argv)
    args.command_line = command_line(argv)
    print(BANNER.format(version=__version__), file=sys.stderr)
    print(f"  Command: {args.command_line}\n", file=sys.stderr)
    t0 = time.time()
    try:
        run(args)
    except QmekuError as exc:
        sys.exit(f"qmeku: {exc}")
    print(f"({time.time() - t0:.0f} s)", file=sys.stderr)
