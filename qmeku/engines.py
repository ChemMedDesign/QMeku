# QMeku: approximate binding energy and per-residue decomposition of docking poses.
# Copyright (C) 2026 SB2CS lab
# SPDX-License-Identifier: GPL-3.0-or-later
# Free software under the GNU General Public License, either version 3 or (at
# your option) any later version; see the LICENSE file.
# Distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
"""Runners for the semiempirical programs: xtb (GFN2-xTB) and MOPAC (PM6-D3H4X).

Each call is a single point without geometry optimisation and returns an :class:`Energy`: the total in kcal/mol
and, when the program prints them, its components (also kcal/mol):

* xtb: ``es`` (isotropic + anisotropic electrostatics), ``xc`` (anisotropic exchange-correlation), ``disp`` (D4
  dispersion), ``rep`` (repulsion), ``gsolv`` (ALPB solvation, 0 in the gas phase) and ``rest`` = total minus
  the rest, the remaining electronic (band) energy. They add up to the total by construction.
* MOPAC: only ``gsolv`` (the COSMO dielectric energy; 0 in the gas phase) and ``rest``. MOPAC does not
  separate electrostatic, dispersion or repulsion terms.

The components of an interaction energy are differences of these terms between the complex and its parts. They
are terms of the GFN2 energy expression, analogous to but not identical with the classical Coulomb and van der
Waals terms; they are not unique and depend on the method.

The programs are found with ``XTB_BIN`` / ``MOPAC_BIN`` or in the PATH (``/opt/mopac/bin/mopac`` as a last
resort).
"""
import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from .ligand import QmekuError

EH_TO_KCAL = 627.509474
EV_TO_KCAL = 23.060548

_XTB_PATTERNS = {
    "total": r":: total energy\s+(-?\d+\.\d+)\s+Eh",
    "es_iso": r":: -> isotropic ES\s+(-?\d+\.\d+)\s+Eh",
    "es_aniso": r":: -> anisotropic ES\s+(-?\d+\.\d+)\s+Eh",
    "xc": r":: -> anisotropic XC\s+(-?\d+\.\d+)\s+Eh",
    "disp": r":: -> dispersion\s+(-?\d+\.\d+)\s+Eh",
    "gsolv": r":: -> Gsolv\s+(-?\d+\.\d+)\s+Eh",
    "rep": r":: repulsion energy\s+(-?\d+\.\d+)\s+Eh",
}


class Energy:
    """Total energy and components, kcal/mol. ``components`` values are floats, or None when the program does not
    print that term."""

    __slots__ = ("total", "components", "etemp")

    def __init__(self, total, components, etemp=None):
        self.total = total
        self.components = components
        self.etemp = etemp      # electronic temperature (K) if the default had to be raised to converge the SCF


def xtb_bin():
    exe = os.environ.get("XTB_BIN") or shutil.which("xtb")
    if not exe:
        raise QmekuError("xtb not found (put it in the PATH or set XTB_BIN)")
    return exe


def mopac_bin():
    exe = os.environ.get("MOPAC_BIN") or shutil.which("mopac") or "/opt/mopac/bin/mopac"
    if not Path(exe).exists():
        raise QmekuError("mopac not found (put it in the PATH or set MOPAC_BIN)")
    return exe


def available():
    """Which programs can be run: {'xtb': bool, 'mopac': bool}."""
    out = {}
    for name, fn in (("xtb", xtb_bin), ("mopac", mopac_bin)):
        try:
            fn()
            out[name] = True
        except QmekuError:
            out[name] = False
    return out


def versions():
    """Version strings of the installed programs, for the log."""
    v = {}
    try:
        p = subprocess.run([xtb_bin(), "--version"], capture_output=True, text=True, timeout=30)
        m = re.search(r"xtb version ([\d.]+)", p.stdout + p.stderr)
        v["xtb"] = m.group(1) if m else "unknown"
    except Exception:
        v["xtb"] = None
    try:
        cff = Path(mopac_bin()).resolve().parents[1] / "CITATION.cff"
        m = re.search(r"^version:\s*(\S+)", cff.read_text(), re.M) if cff.exists() else None
        if not m:
            p = subprocess.run([mopac_bin(), "--version"], capture_output=True, text=True, timeout=15)
            m = re.search(r"version\s+([\d.]+)", p.stdout + p.stderr, re.I)
        v["mopac"] = m.group(1) if m else "unknown"
    except Exception:
        v["mopac"] = None
    return v


def _check_timeout(timeout):
    if timeout is not None and timeout <= 0:
        raise QmekuError(f"timeout must be > 0 (got {timeout})")


def run_xtb(atoms, charge, solvent="water", gfn=2, threads=1, keep_dir=None, timeout=None):
    """GFN-xTB single point of ``atoms`` = [(element, (x, y, z), label)]; ``solvent=None`` is the gas phase."""
    _check_timeout(timeout)
    deadline = None if timeout is None else time.monotonic() + timeout
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(keep_dir) if keep_dir else Path(tmp)
        d.mkdir(parents=True, exist_ok=True)
        lines = [str(len(atoms)), "qmeku"] + [f"{el} {x:.6f} {y:.6f} {z:.6f}" for el, (x, y, z), _ in atoms]
        (d / "in.xyz").write_text("\n".join(lines) + "\n")
        cmd = [xtb_bin(), "in.xyz", "--gfn", str(gfn), "--chrg", str(charge), "--sp"]
        if solvent:
            cmd += ["--alpb", solvent]
        # xtb segfaults on systems of several hundred atoms with the default 8 MB stack: raise the limit in a
        # shell wrapper (preexec_fn is not safe with threads)
        cmd = ["bash", "-c", 'ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -H -s)" 2>/dev/null; exec "$@"',
               "xtb-wrapper"] + cmd
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads), OMP_STACKSIZE="1G")
        # The SCF of some systems (e.g. an isolated ion pair in the gas phase) does not converge with the default
        # settings; retry with more cycles and a higher electronic temperature (a few tenths of kcal/mol at most on
        # systems that converge both ways) and record it in the result.
        etemp_used = None
        for extra, etemp in (([], None), (["--cycles", "2000", "--etemp", "600"], 600),
                             (["--cycles", "2000", "--etemp", "1000"], 1000), (["--cycles", "2000", "--etemp", "2000"], 2000)):
            for stale in ("xtbrestart", "charges"):
                (d / stale).unlink(missing_ok=True)
            rem = None
            if deadline is not None:
                rem = deadline - time.monotonic()
                if rem <= 0:
                    raise QmekuError(f"xtb timed out after {timeout} s")
            try:
                proc = subprocess.run(cmd + extra, cwd=d, env=env, capture_output=True, text=True, timeout=rem)
            except subprocess.TimeoutExpired:
                raise QmekuError(f"xtb timed out after {timeout} s") from None
            text = proc.stdout + proc.stderr
            if proc.returncode == 0 and "normal termination" in text:
                etemp_used = etemp
                break
            if "did not converge" not in text:
                raise QmekuError("xtb failed: " + text[-500:])
        else:
            raise QmekuError("xtb: the SCF did not converge even with an electronic temperature of 2000 K: " + text[-300:])
        val = {}
        for key, pat in _XTB_PATTERNS.items():
            m = re.search(pat, proc.stdout)
            val[key] = float(m.group(1)) * EH_TO_KCAL if m else None
        if val["total"] is None:
            raise QmekuError("xtb: no total energy found: " + text[-300:])
        comp = {"es": (val["es_iso"] or 0.0) + (val["es_aniso"] or 0.0), "xc": val["xc"] or 0.0,
                "disp": val["disp"] or 0.0, "rep": val["rep"] or 0.0, "gsolv": val["gsolv"] or 0.0}
        comp["rest"] = val["total"] - sum(comp.values())
        return Energy(val["total"], comp, etemp_used)


def run_mopac(atoms, charge, method="PM6-D3H4X", eps=78.4, mozyme=False, threads=1, keep_dir=None, timeout=None):
    """MOPAC single point (1SCF) of ``atoms``; the energy is the heat of formation. ``eps=None`` is the gas phase."""
    _check_timeout(timeout)
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(keep_dir) if keep_dir else Path(tmp)
        d.mkdir(parents=True, exist_ok=True)
        kw = [method, "1SCF", f"CHARGE={charge}", f"THREADS={threads}", "GEO-OK"]
        if eps:
            kw.append(f"EPS={eps}")
        if mozyme:
            kw.append("MOZYME")
        lines = [" ".join(kw), "qmeku", ""]
        lines += [f"{el} {x:.6f} 0 {y:.6f} 0 {z:.6f} 0" for el, (x, y, z), _ in atoms]
        (d / "in.mop").write_text("\n".join(lines) + "\n")
        try:
            proc = subprocess.run([mopac_bin(), "in.mop"], cwd=d, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise QmekuError(f"mopac timed out after {timeout} s") from None
        out = (d / "in.out").read_text() if (d / "in.out").exists() else ""
        m = re.search(r"FINAL HEAT OF FORMATION\s*=\s*(-?\d+\.\d+)\s+KCAL/MOL", out)
        if not m or "JOB ENDED NORMALLY" not in out:
            raise QmekuError("mopac failed: " + (out[-600:] or proc.stdout + proc.stderr))
        total = float(m.group(1))
        md = re.search(r"DIELECTRIC ENERGY\s*=\s*(-?\d+\.\d+)\s+EV", out)
        gsolv = float(md.group(1)) * EV_TO_KCAL if md else 0.0
        return Energy(total, {"es": None, "xc": None, "disp": None, "rep": None, "gsolv": gsolv, "rest": total - gsolv})


def make_runner(method, solvent_on, solvent="water", timeout=None):
    """Function (atoms, charge) -> Energy for ``method`` ('xtb' or 'mopac') with or without implicit solvent."""
    if method == "xtb":
        return lambda atoms, q: run_xtb(atoms, q, solvent=solvent if solvent_on else None, threads=1, timeout=timeout)
    if method == "mopac":
        if solvent_on and solvent != "water":
            raise QmekuError(f"MOPAC solvent must be 'water' (got {solvent!r})")
        return lambda atoms, q: run_mopac(atoms, q, eps=78.4 if solvent_on else None, threads=1, timeout=timeout)
    raise QmekuError(f"unknown method {method!r} (use xtb or mopac)")
