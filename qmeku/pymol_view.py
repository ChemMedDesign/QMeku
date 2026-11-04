# QMeku: approximate binding energy and per-residue decomposition of docking poses.
# Copyright (C) 2026 SB2CS lab
# SPDX-License-Identifier: GPL-3.0-or-later
# Free software under the GNU General Public License, either version 3 or (at
# your option) any later version; see the LICENSE file.
# Distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
"""PyMOL view of the per-residue decomposition (written by default next to the CSV files; ``--no-pymol`` skips it).

For the best poses (lowest E_bind of the first method) and for every method it writes, in ``<outdir>/pymol/``:

* ``<pose>_<method>.pml``  self-contained PyMOL script (plain commands, relative paths) that rebuilds the scene;
* ``<pose>_<method>.png``  image, only if PyMOL is found (``pymol`` in the PATH or the
  ``PYMOL_BIN`` variable); otherwise the ``.pml`` is still written and can be opened with ``pymol file.pml``;
* ``receptor.pdb``, ``ligand_<pose>.sdf``  the coordinates the scripts load, and
  ``<pose>_<method>_receptor_Bfactor.pdb`` (energy of each residue in the B-factor column, with PyMOL).

Scene: every residue in contact is coloured by its contribution ``{method}_pair`` (kcal/mol): most negative = red,
close to 0 = light grey, positive = light blue (only if there are unfavourable residues). The ligand is drawn in
pale yellow (small spheres + sticks) and only the polar hydrogens of the pocket and the ligand are shown. The most
important residues are labelled (``Asp189``) and the polar contacts of those residues with the ligand are dashed.
With several poses an extra scene uses the mean contribution over poses (``residues_summary.csv``).
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

GREY = (0.85, 0.85, 0.85)          # contribution ~ 0
RED = (0.90, 0.00, 0.00)           # most negative contribution
LIGHT_BLUE = (0.55, 0.78, 1.00)    # positive (unfavourable) contribution


def _mix(c0, c1, t):
    t = max(0.0, min(1.0, t))
    return tuple(a + (b - a) * t for a, b in zip(c0, c1))


def energy_color(e, emin, emax, gamma=0.6):
    """RGB of a contribution ``e``: ``emin`` (<0) red, 0 light grey, ``emax`` (>0) light blue. ``gamma`` < 1
    emphasises small contributions next to a dominant one (1 = linear)."""
    if e < 0 and emin < 0:
        return _mix(GREY, RED, (e / emin) ** gamma)
    if e > 0 and emax > 0:
        return _mix(GREY, LIGHT_BLUE, (e / emax) ** gamma)
    return GREY


def parse_label(label):
    """'ASP189', 'A:ASP189', 'ALA221A' -> (chain or None, resname, resi with insertion code)."""
    chain = None
    if ":" in label:
        chain, label = label.split(":", 1)
    m = re.match(r"^([A-Za-z0-9]{3})(-?\d+)([A-Za-z]?)$", label)
    if not m:
        raise ValueError("unrecognised residue label: " + label)
    return chain, m.group(1), m.group(2) + m.group(3)


def _selection(r):
    s = "receptor and polymer and resi %s" % r["resi"]
    return "(" + s + (" and chain %s" % r["chain"] if r["chain"] else "") + ")"


def _rows_from_pose(res, method, mol, contact_max):
    """Residue entries of one pose result (energy, contact flag, polar contact with the ligand)."""
    out = []
    for w in res["residues"]:
        chain, resname, resi = parse_label(w["residue"])
        row = {"chain": chain, "resname": resname, "resi": resi, "e": w[method]["pair"], "contact": bool(w["contact"]),
               "lig_idx": None, "rec_atom": None}
        li, ra, d = w.get("lig_heavy_atom"), w.get("rec_heavy_atom"), w.get("min_dist_heavy")
        if li is not None and ra and d is not None and d <= contact_max:
            if mol.GetAtomWithIdx(int(li) - 1).GetSymbol() in ("N", "O", "S") and ra[:1] in ("N", "O", "S"):
                row["lig_idx"], row["rec_atom"] = int(li), ra           # qmeku's atom index is 1-based = PyMOL index
        out.append(row)
    return out


def _rows_from_summary(summ, method):
    out = []
    for s in summ:
        chain, resname, resi = parse_label(s["residue"])
        out.append({"chain": chain, "resname": resname, "resi": resi, "e": s[method + "_mean"],
                    "contact": s["n_contact"] > 0, "lig_idx": None, "rec_atom": None})
    return out


def scene_commands(rows, ligand_file, emin, emax, gamma=0.6, top=8):
    """Plain PyMOL commands of one scene (no Python, relative paths)."""
    L = ["reinitialize", "load receptor.pdb, receptor", "load %s, ligand" % ligand_file,
         "hide everything", "bg_color white",
         "show cartoon, receptor and polymer", "color grey90, receptor",
         "set cartoon_transparency, 0.7, receptor", "set cartoon_fancy_helices, 1",
         "show sticks, ligand", "show spheres, ligand", "color paleyellow, ligand", "util.cnc ligand",
         "set sphere_scale, 0.2, ligand", "set sphere_scale, 0.1, ligand and hydro",
         "set stick_radius, 0.14, receptor", "set stick_radius, 0.18, ligand",
         "set cartoon_loop_radius, 0.12", "set cartoon_oval_length, 0.9", "set cartoon_rect_length, 0.9"]
    active = [r for r in rows if r["contact"]]
    for k, r in enumerate(active):
        name = "qm_%03d" % k
        sel = _selection(r)
        L += ["set_color %s, [%.3f, %.3f, %.3f]" % ((name,) + energy_color(r["e"], emin, emax, gamma)),
              "show sticks, %s" % sel, "util.cnc %s" % sel, "color %s, %s and elem C" % (name, sel),
              "set cartoon_transparency, 0.0, %s" % sel]
    L.append("alter receptor, b=0.0")                       # energy per residue in the B factor
    L += ["alter %s, b=%.4f" % (_selection(r), r["e"]) for r in rows]
    L.append("alter ligand, b=0.0")
    L += ["set label_font_id, 7", "set label_size, 13", "set label_color, black",
          "set label_outline_color, white", "set label_position, [1.5, 1.5, 4]"]
    n = 0
    for r in sorted((x for x in rows if x["e"] != 0), key=lambda x: abs(x["e"]), reverse=True)[:top]:
        L.append('label (%s and name CA), "%s"' % (_selection(r), r["resname"].capitalize() + r["resi"]))
        if r["lig_idx"] is not None:
            n += 1
            L.append("distance qm_polar_%d, ligand and index %d, %s and name %s" %
                     (n, r["lig_idx"], _selection(r), r["rec_atom"]))
    L += ["set dash_gap, 0.25", "set dash_radius, 0.07", "set dash_color, black", "set dash_round_ends, 0"]
    if n:
        L += ["hide labels, qm_polar_*"]
    L += ["select qm_pocket, " + " or ".join(["ligand"] + [_selection(r) for r in active]),
          "orient qm_pocket", "zoom qm_pocket, 4",
          "set ray_shadows, 0", "set antialias, 2", "set depth_cue, 0", "set specular, 0.2",
          "hide everything, (receptor or ligand) and hydro and not neighbor elem N+O+S", "deselect"]
    return L


def find_pymol():
    return os.environ.get("PYMOL_BIN") or shutil.which("pymol")


def _scale(rows_list, vmin=None, vmax=None):
    es = [r["e"] for rows in rows_list for r in rows]
    neg, pos = [e for e in es if e < 0], [e for e in es if e > 0]
    emin = vmin if vmin is not None else (min(neg) if neg else -1.0)
    emax = vmax if vmax is not None else (max(pos) if pos else 0.0)
    return (emin if emin < 0 else -1.0), emax


def _safe(name):
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", str(name))


def write_scenes(outdir, receptor_text, results, methods, pose_mols, summary=None, n_poses=3, top=8,
                 contact_max=3.5, gamma=0.6, make_images=True, timeout=900, lang="en"):
    """Write the PyMOL files in ``<outdir>/pymol``. ``results`` = [(pose_name, result)], ``pose_mols`` =
    {pose_name: (rdkit_mol, conformer_id)}, ``summary`` = rows of :func:`qmeku.report.summarize_residues`.
    Returns a list of messages for the user."""
    from rdkit import Chem
    pdir = Path(outdir) / "pymol"
    first = methods[0]
    ranked = sorted((x for x in results if x[1].get("residues")), key=lambda x: x[1]["methods"][first]["e_bind"])
    chosen = ranked if n_poses == 0 else ranked[:n_poses]
    if not chosen:
        return []
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "receptor.pdb").write_text(receptor_text)

    sdf = {}

    def ligand_file(pose):
        if pose not in sdf:
            mol, cid = pose_mols[pose]
            fn = "ligand_%s.sdf" % _safe(pose)
            (pdir / fn).write_text(Chem.MolToMolBlock(mol, confId=cid) + "$$$$\n")
            sdf[pose] = fn
        return sdf[pose]

    scenes = []                                              # (tag, commands, description)
    for m in methods:
        data = [(p, _rows_from_pose(r, m, pose_mols[p][0], contact_max)) for p, r in chosen]
        emin, emax = _scale([rows for _, rows in data])      # same scale for all poses of a method
        for p, rows in data:
            scenes.append(("%s_%s" % (_safe(p), m),
                           scene_commands(rows, ligand_file(p), emin, emax, gamma, top),
                           "pose %s, %s_pair (kcal/mol), scale %.2f .. 0%s" % (
                               p, m, emin, " .. +%.2f" % emax if emax > 0 else "")))
        if summary:
            rows = _rows_from_summary(summary, m)
            e0, e1 = _scale([rows])
            scenes.append(("summary_%s" % m, scene_commands(rows, ligand_file(chosen[0][0]), e0, e1, gamma, top),
                           "mean over poses, %s_mean (kcal/mol); ligand = %s" % (m, chosen[0][0])))

    for tag, cmds, info in scenes:
        (pdir / (tag + ".pml")).write_text(
            "# QMeku -> PyMOL | %s\n# red = most negative, light grey ~ 0, light blue = positive\n"
            "# Run from this folder:  pymol %s.pml\n" % (info, tag) + "\n".join(cmds) + "\n")

    msgs = [("Archivos de PyMOL (%d escenas .pml) en %s/" if lang == "es" else "PyMOL files (%d .pml scenes) in %s/")
            % (len(scenes), pdir)]
    exe = find_pymol() if make_images else None
    if exe:
        run = pdir / "_render.pml"
        run.write_text("\n".join("@%s.pml\nsave %s_receptor_Bfactor.pdb, receptor\n"
                                 "png %s.png, 2000, 1500, dpi=300, ray=1" % (t, t, t) for t, _, _ in scenes)
                       + "\nquit\n")
        try:
            subprocess.run([exe, "-cq", run.name], cwd=pdir, capture_output=True, timeout=timeout)
            msgs.append("PyMOL: .png" if all((pdir / (t + ".png")).exists() for t, _, _ in scenes)
                        else "PyMOL: some .png were not written (open the .pml files)")
        except (subprocess.TimeoutExpired, OSError) as exc:
            msgs.append("PyMOL: %s" % str(exc)[:100])
        finally:
            run.unlink(missing_ok=True)
    elif make_images:
        msgs.append("PyMOL no encontrado (PATH o PYMOL_BIN): solo se escribieron los .pml; ábrelos con 'pymol archivo.pml'"
                    if lang == "es" else
                    "PyMOL not found (PATH or PYMOL_BIN): only the .pml files were written; open with 'pymol file.pml'")
    return msgs
