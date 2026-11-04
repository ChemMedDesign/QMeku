"""Binding-energy approximation and per-residue decomposition of docking poses with xtb and/or MOPAC.

For one pose (a ligand in a receptor with hydrogens) it builds a pocket model (residues within ``cutoff`` of the
ligand, cut peptide bonds capped with hydrogen link atoms, see :mod:`qmeku.structure`) and runs single points
without geometry optimisation, in the implicit solvent and (by default) in the gas phase:

* ``e_bind``  E(pocket+ligand) - E(pocket) - E(ligand), in the implicit solvent: the approximate binding energy
              of the pose. It does not include the entropy or the energy of reorganising the ligand or the protein.
* ``e_gas``   the same in the gas phase, and ``d_solv = e_bind - e_gas`` the change of the solvation energy
              upon binding (positive = desolvation penalty).
* per residue ``pair`` = E(ligand+F) - E(ligand) - E(F), F being the residue capped on both sides: the interaction
              of the ligand with that residue alone (also ``pair_gas`` and ``dsolv`` = ``pair - pair_gas``).
* components (xtb): the same differences of the electrostatic, dispersion, repulsion, XC, solvation and
              remaining-electronic terms of the GFN2 energy (see :mod:`qmeku.engines`).

The pair energies do not add up to ``e_bind`` and the balance is reported explicitly: in the gas phase because
of the cooperative (many-body) effects, and in the solvent because each pair is solvated on its own.
"""
import math
import os
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from . import structure as pk
from .engines import make_runner
from .ligand import QmekuError, build_complex

STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "CYM", "GLN", "GLU", "GLY", "HID", "HIE", "HIP", "ILE", "LEU",
            "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH", "LYN", "ACE", "NME"}
TERMS = ("es", "xc", "disp", "rep", "gsolv", "rest")


class InputError(QmekuError):
    """The receptor or the ligand cannot be used (protonation, unsupported residues, charge)."""


def check_residues(cx, keys):
    """The pocket must contain only amino acids with an explicit protonation state and valid termini/disulfides."""
    for k in keys:
        r = cx.residues[k]
        name = r["resname"]
        if name == "HIS":
            raise InputError(f"{cx.label(k)}: generic HIS has no protonation state; use HID/HIE/HIP "
                             "(prepare the receptor with PDB2PQR)")
        if name not in STANDARD:
            raise InputError(f"{cx.label(k)}: residue {name} is not supported in the pocket (only amino acids; "
                             "remove waters, ions and cofactors from the receptor)")
        if name == "CYX" and k not in cx.disulfides:
            raise InputError(f"{cx.label(k)}: CYX has no disulfide partner within 2.5 A")
        if "OXT" in r["atoms"] and pk._neighbour(cx, k, +1) is not None:
            raise InputError(f"{cx.label(k)}: internal peptide-bonded residue contains a terminal OXT atom")
        if pk._neighbour(cx, k, -1) is not None and any(h in r["atoms"] for h in ("H1", "H2", "H3")):
            raise InputError(f"{cx.label(k)}: internal peptide-bonded residue contains N-terminal H1/H2/H3 atoms")


def closest_contact(res_atoms, lig_atoms, heavy_only=False):
    """(distance, receptor atom name, ligand atom index (1-based), ligand atom name) of the closest pair of atoms,
    hydrogens included or (``heavy_only``) only the heavy atoms."""
    best, d_best = (float("nan"), "", 0, ""), 1e9
    for rn, (el, xyz) in res_atoms.items():
        if heavy_only and el == "H":
            continue
        for i, (lel, lxyz, label) in enumerate(lig_atoms):
            if heavy_only and lel == "H":
                continue
            d = math.dist(xyz, lxyz)
            if d < d_best:
                d_best, best = d, (d, rn, i + 1, label.split(":", 1)[1])
    return best


def _diff(e_cx, e_a, e_b):
    """Interaction energy and its components: complex - part a - part b (Energy objects)."""
    total = e_cx.total - e_a.total - e_b.total
    comp = {}
    for t in TERMS:
        vals = (e_cx.components.get(t), e_a.components.get(t), e_b.components.get(t))
        comp[t] = None if any(v is None for v in vals) else vals[0] - vals[1] - vals[2]
    return total, comp


def analyze_pose(complex_pdb, ligand_charge, methods=("xtb",), solvent="water", gas=True, cutoff=4.0, pairs=True,
                 workers=None, ligand_resname="LIG", residue_keys=None, concap=True, timeout=None):
    """Analyse one pose given as a complex PDB path (receptor ATOM records + the ligand as HETATM residue).

    Returns ``{"pocket", "ligand", "seconds", "residues": [...], "methods": {name: {...}}}``. Each residue row has the
    geometry (``min_dist`` and the two atoms, hydrogens included, and ``min_dist_heavy`` with the heavy atoms)
    and, for each method, ``row[method] = {pair, pair_gas, dsolv, comp, comp_gas}``. ``residue_keys`` fixes the
    pocket (residue keys of :class:`qmeku.structure.Complex`); by default, the residues within ``cutoff`` of the
    ligand in this pose."""
    cx = pk.Complex(str(complex_pdb), ligand_resname)
    if residue_keys is None:
        sel = pk.select_residues(cx, cutoff)
    else:
        wanted = set(residue_keys)
        sel = [k for k in cx.protein_keys if k in wanted]
    if not sel:
        raise InputError(f"no receptor residue within {cutoff} A of the ligand")
    P, used = pk.pocket_atoms(cx, sel)
    check_residues(cx, used)
    L = cx.ligand_atoms()
    qP = pk.system_charge(cx, used)
    try:
        for k in used:
            qF_chk = pk.residue_charge(cx, k)
            if pk.n_electrons(pk.fragment_atoms(cx, k), qF_chk) % 2:
                raise InputError(f"odd number of electrons in residue {cx.label(k)} (charge {qF_chk}): "
                                 "check the protonation and termini of the receptor")
        for name, atoms, q in (("pocket", P, qP), ("pocket + ligand", P + L, qP + ligand_charge), ("ligand", L, ligand_charge)):
            if pk.n_electrons(atoms, q) % 2:
                raise InputError(f"odd number of electrons in the {name} (charge {q}): check the protonation of the "
                                 "receptor and the formal charge of the ligand")
    except ValueError as exc:
        raise InputError(str(exc)) from None
    min_d = pk.link_atoms_distance_check(P + L)
    if min_d < 0.5:
        raise InputError(f"overlapping atoms in the pocket/ligand model (minimum distance {min_d:.3f} A < 0.5 A)")

    solv_on = solvent not in (None, "none", "gas")
    gas = bool(gas and solv_on)          # without an implicit solvent the single phase is the gas phase itself
    phases = [("solv", solv_on)] + ([("gas", False)] if gas else [])
    runners = {(m, ph): make_runner(m, on, solvent, timeout=timeout) for m in methods for ph, on in phases}
    systems = {("L",): (L, ligand_charge), ("P",): (P, qP), ("PL",): (P + L, qP + ligand_charge)}
    if pairs:
        for k in sel:
            F = pk.fragment_atoms(cx, k)
            qF = pk.residue_charge(cx, k)
            systems[("F", k)] = (F, qF)
            systems[("FL", k)] = (F + L, qF + ligand_charge)
    bonds = {}
    if pairs and concap:
        for k, k2 in pk.peptide_bonds(cx, sel):
            atoms = pk.concap_atoms(cx, k, k2)
            if atoms is not None:
                bonds[(k, k2)] = atoms
                systems[("C", k, k2)] = (atoms, 0)
                systems[("CL", k, k2)] = (atoms + L, ligand_charge)
    jobs = [((m, ph) + key, runners[(m, ph)], atoms, q) for (m, ph) in runners for key, (atoms, q) in systems.items()]
    t0 = time.time()
    pool = ThreadPoolExecutor(max_workers=workers or os.cpu_count() or 1)
    try:
        futures = [pool.submit(lambda j=j: (j[0], j[1](j[2], j[3]))) for j in jobs]
        energies = {}
        for fut in as_completed(futures):
            k, val = fut.result()
            energies[k] = val
    except Exception:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    else:
        pool.shutdown(wait=True)
    fallbacks = [(j[0], energies[j[0]].etemp) for j in jobs if energies[j[0]].etemp]

    geometry = {}
    if pairs:
        for k in sel:
            r = cx.residues[k]
            d, ra, li, ln = closest_contact(r["atoms"], L)
            dh, rh, lhi, lhn = closest_contact(r["atoms"], L, heavy_only=True)
            near = sum(1 for _, (el, xyz) in r["atoms"].items() if min(math.dist(xyz, l[1]) for l in L) <= cutoff)
            geometry[k] = {"residue": cx.label(k), "resname": r["resname"], "resnum": r["resnum"], "index": cx.index(k),
                           "min_dist": d, "lig_atom": li, "lig_atom_name": ln, "rec_atom": ra,
                           "min_dist_heavy": dh, "lig_heavy_atom": lhi, "lig_heavy_name": lhn, "rec_heavy_atom": rh,
                           "n_atoms_close": near, "contact": d <= cutoff}
    out_methods = {}
    for m in methods:
        res = {}
        for ph, _ in phases:
            def e(*key, ph=ph):
                return energies[(m, ph) + key]
            tot, comp = _diff(e("PL"), e("P"), e("L"))
            res[ph] = {"e": tot, "comp": comp}
            if pairs:
                raw = {k: _diff(e("FL", k), e("L"), e("F", k)) for k in sel}
                cc = {b: _diff(e("CL", *b), e("L"), e("C", *b)) for b in bonds}
                adj = {}
                for k in sel:
                    t, c = raw[k]
                    c = dict(c)
                    for b, (ct, cco) in cc.items():
                        if k in b:
                            t -= 0.5 * ct
                            for term in TERMS:
                                if c[term] is not None and cco[term] is not None:
                                    c[term] -= 0.5 * cco[term]
                    adj[k] = (t, c)
                res[ph]["pairs"] = adj
                res[ph]["pairs_raw"] = {k: raw[k][0] for k in sel}
                res[ph]["concap_sum"] = sum(v[0] for v in cc.values())
        e_bind = res["solv"]["e"]
        e_gas = res["gas"]["e"] if gas else (e_bind if not solv_on else None)
        entry = {"e_bind": e_bind, "e_gas": e_gas, "d_solv": e_bind - res["gas"]["e"] if gas else None,
                 "comp": res["solv"]["comp"] if solv_on else None,
                 "comp_gas": res["gas"]["comp"] if gas else (res["solv"]["comp"] if not solv_on else None)}
        if pairs:
            ps = {k: res["solv"]["pairs"][k][0] for k in sel}
            pg = {k: res["gas"]["pairs"][k][0] for k in sel} if gas else (ps if not solv_on else None)
            entry["pair_sum"] = sum(ps.values())
            entry["concap_sum"] = res["solv"]["concap_sum"]
            entry["pair_sum_raw"] = sum(res["solv"]["pairs_raw"].values())
            entry["residual"] = e_bind - entry["pair_sum"]
            if pg is not None:
                entry["pair_sum_gas"] = sum(pg.values())
                entry["residual_gas"] = e_gas - entry["pair_sum_gas"]
            if gas:
                entry["pair_dsolv_sum"] = entry["pair_sum"] - entry["pair_sum_gas"]
                entry["residual_solv"] = entry["d_solv"] - entry["pair_dsolv_sum"]
            entry["per_residue"] = {
                k: {"pair": ps[k], "pair_raw": res["solv"]["pairs_raw"][k],
                    "pair_gas": pg[k] if pg is not None else None,
                    "dsolv": ps[k] - pg[k] if gas else None,
                    "comp": res["solv"]["pairs"][k][1] if solv_on else None,
                    "comp_gas": res["gas"]["pairs"][k][1] if gas else (res["solv"]["pairs"][k][1] if not solv_on else None)}
                for k in sel}
        out_methods[m] = entry
    rows = []
    if pairs:
        for k in sel:
            row = dict(geometry[k])
            for m in methods:
                row[m] = out_methods[m]["per_residue"][k]
            rows.append(row)
    return {"pocket": {"selected": len(sel), "residues": len(used), "atoms": len(P), "charge": qP},
            "ligand": {"atoms": len(L), "charge": ligand_charge}, "seconds": time.time() - t0,
            "solvent": solvent if solv_on else "none", "solvent_on": solv_on, "gas": gas, "cutoff": cutoff, "scf_fallbacks": [(" ".join(str(x if isinstance(x, (str, int)) else cx.label(x)) for x in key), t) for key, t in fallbacks],
            "concap": bool(bonds) if pairs else False, "n_concaps": len(bonds), "residues": rows,
            "methods": out_methods}


def pocket_residues(receptor, mol, cutoff=4.0, conf_ids=None):
    """Residue keys within ``cutoff`` of the ligand in any of its conformers (the union is the pocket shared by
    several poses, needed to compare their energies)."""
    keys = []
    for cid in (conf_ids if conf_ids is not None else [c.GetId() for c in mol.GetConformers()]):
        text = build_complex(receptor, mol, [cid])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "complex.pdb"
            path.write_text(text)
            cx = pk.Complex(str(path), "LIG")
            keys += [k for k in pk.select_residues(cx, cutoff) if k not in keys]
    return keys


def analyze_ligand(receptor, mol, ligand_charge=None, name="pose", conf_ids=None, timeout=None, **kw):
    """Analyse an RDKit ligand (hydrogens and 3D coordinates; one pose per conformer) against a receptor PDB
    (path or text). Returns ``[(pose_name, result)]``. ``kw`` are the arguments of :func:`analyze_pose`."""
    from rdkit import Chem
    if ligand_charge is None:
        ligand_charge = Chem.GetFormalCharge(mol)
    confs = list(conf_ids) if conf_ids is not None else [c.GetId() for c in mol.GetConformers()]
    out = []
    for cid in confs:
        text = build_complex(receptor, mol, [cid])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "complex.pdb"
            path.write_text(text)
            res = analyze_pose(path, ligand_charge, timeout=timeout, **kw)
        out.append((name if len(confs) == 1 else f"{name}_c{cid}", res))
    return out
