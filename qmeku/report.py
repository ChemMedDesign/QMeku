"""Output files and log of QMeku.

Files (``write_outputs``): ``poses.csv``, ``residues.csv``, ``residues_summary.csv`` (several poses) and
``qmeku.log``. The log (``format_log``) has the binding energy of every pose, the balance between gas phase and
solvation, the GFN2 energy terms (electrostatics, dispersion, repulsion, XC, solvation), and a per-residue table
that reconciles with the binding energy.
"""
import csv
import datetime
import statistics
from pathlib import Path

TERMS = ("es", "xc", "disp", "rep", "gsolv", "rest")
LABEL = {"xtb": "GFN2-xTB", "mopac": "PM6-D3H4X"}

TEXT = {
    "es": {
        "title": "QMeku {ver}: energía de unión aproximada y descomposición por residuo",
        "date": "Fecha", "receptor": "Receptor", "poses": "Poses", "methods": "Métodos",
        "pocket": "Bolsillo", "pocket_common": "residuos a menos de {c} Å del ligando en alguna pose (los mismos para todas)",
        "pocket_each": "residuos a menos de {c} Å del ligando en cada pose",
        "warn": ("AVISO: puntos únicos sin optimizar la geometría. E_bind es la energía de interacción ligando–bolsillo en "
                 "solvente implícito; no es una energía libre (no incluye entropía ni la reorganización del ligando o de la proteína) "
                 "y no está validada frente a datos experimentales."),
        "pose": "Pose", "ligand_charge": "carga del ligando", "pocket_line": "bolsillo {r} residuos / {a} átomos, carga {q:+d}",
        "bind": "ENERGÍA DE UNIÓN APROXIMADA  E_bind = E(complejo) − E(bolsillo) − E(ligando), en solvente implícito (kcal/mol)",
        "bind_gas": "ENERGÍA DE INTERACCIÓN EN FASE GAS  E_bind = E(complejo) − E(bolsillo) − E(ligando), fase gas (kcal/mol)",
        "balance": "Balance de E_bind", "gas": "E_int en fase gas",
        "dsolv": "cambio de solvatación (penalización por desolvatación)", "ebind": "E_bind (solvente)",
        "terms_gas": "Términos de E_int en fase gas (energía de interacción, kcal/mol)",
        "t_es": "electrostática (ES)", "t_disp": "dispersión D4", "t_rep": "repulsión", "t_xc": "XC anisotrópico",
        "t_rest": "resto electrónico", "t_gsolv": "solvatación (Gsolv)", "t_total": "total",
        "terms": "Términos de GFN2-xTB en solvente implícito (energía de interacción, kcal/mol)",
        "terms_note": ("Son términos de la expresión de la energía de GFN2: electrostática (ES), dispersión D4, "
                       "repulsión, XC anisotrópico, solvatación ALPB (Gsolv) y resto electrónico; su suma es E_bind. "
                       "No son particiones termodinámicas únicas ni equivalen a términos de campo de fuerza clásico. "
                       "MOPAC no separa estos términos."),
        "perres": "Contribución por residuo: energía de interacción del ligando con cada residuo por separado (kcal/mol)",
        "res": "residuo", "dist": "dist.", "atoms": "átomos", "heavy": "pesados",
        "sum": "Σ residuos", "rest": "resto no aditivo",
        "more": "(+{n} residuos más en residues.csv)",
        "split": "Desglose por residuo: fase gas y cambio de solvatación (kcal/mol); términos de GFN2-xTB en fase gas",
        "why": ("El resto no aditivo reúne dos efectos de signo contrario: los efectos cooperativos entre residuos (fase gas) y la "
                "solvatación, porque cada par se calcula con su propio disolvente implícito e incluye una penalización de desolvatación "
                "que por separado se sobrestima. No es solo la desolvatación."),
        "concap": ("Corrección de concaps: {n} enlaces peptídicos entre residuos del bolsillo; el grupo peptídico compartido por dos "
                   "fragmentos vecinos se resta una vez para no contarlo dos veces. Σ residuos sin corregir → corregida: {txt}"),
        "fallback": ("AVISO: {n} cálculo(s) de xtb no convergieron con los parámetros por defecto y se repitieron con temperatura "
                     "electrónica mayor (efecto típico de décimas de kcal/mol): {txt}"),
        "concap_off": "Sin corrección de concaps (--no-concap): los grupos peptídicos compartidos se cuentan en dos residuos.",
        "summary": "RESUMEN DE {n} POSES", "rank": "ordenadas por E_bind de {m}", "pose_col": "pose",
        "res_mean": "Contribución media por residuo ({m}, kcal/mol)", "contact": "contacto", "mean": "media",
        "median": "mediana", "sd": "desv.", "best": "mejor dist.",
        "filtered": "(las {k} poses con E_bind ≤ {x:.1f}; las demás se excluyen del promedio por residuo)",
        "files": "Archivos escritos en", "failed": "Poses que fallaron",
        "noncommon": "Aviso: bolsillos distintos por pose; las energías de poses distintas no son estrictamente comparables.",
    },
    "en": {
        "title": "QMeku {ver}: approximate binding energy and per-residue decomposition",
        "date": "Date", "receptor": "Receptor", "poses": "Poses", "methods": "Methods",
        "pocket": "Pocket", "pocket_common": "residues within {c} A of the ligand in any pose (the same for all)",
        "pocket_each": "residues within {c} A of the ligand in each pose",
        "warn": ("WARNING: single points without geometry optimisation. E_bind is the ligand-pocket interaction energy in "
                 "implicit solvent; it is not a free energy (no entropy, no reorganisation of the ligand or the protein) "
                 "and it has not been validated against experiment."),
        "pose": "Pose", "ligand_charge": "ligand charge", "pocket_line": "pocket {r} residues / {a} atoms, charge {q:+d}",
        "bind": "APPROXIMATE BINDING ENERGY  E_bind = E(complex) - E(pocket) - E(ligand), implicit solvent (kcal/mol)",
        "bind_gas": "GAS-PHASE INTERACTION ENERGY  E_bind = E(complex) - E(pocket) - E(ligand), gas phase (kcal/mol)",
        "balance": "Balance of E_bind", "gas": "E_int in the gas phase",
        "dsolv": "solvation change (desolvation penalty)", "ebind": "E_bind (solvent)",
        "terms_gas": "Terms of the gas-phase E_int (interaction energy, kcal/mol)",
        "t_es": "electrostatics (ES)", "t_disp": "D4 dispersion", "t_rep": "repulsion", "t_xc": "anisotropic XC",
        "t_rest": "electronic remainder", "t_gsolv": "solvation (Gsolv)", "t_total": "total",
        "terms": "GFN2-xTB terms in implicit solvent (interaction energy, kcal/mol)",
        "terms_note": ("These are terms of the GFN2 energy expression: electrostatics (ES), D4 dispersion, repulsion, "
                       "anisotropic XC, ALPB solvation (Gsolv) and electronic remainder; they add up to E_bind. "
                       "They are not unique thermodynamic partitions and are not equivalent to classical force-field "
                       "terms. MOPAC does not separate these terms."),
        "perres": "Per-residue contribution: interaction energy of the ligand with each residue alone (kcal/mol)",
        "res": "residue", "dist": "dist.", "atoms": "atoms", "heavy": "heavy",
        "sum": "sum of residues", "rest": "non-additive remainder",
        "more": "(+{n} more residues in residues.csv)",
        "split": "Per-residue split: gas phase and solvation change (kcal/mol); GFN2-xTB terms in the gas phase",
        "why": ("The non-additive remainder joins two effects of opposite sign: the cooperative effects between residues (gas phase) and "
                "the solvation, because each pair is calculated with its own implicit solvent and carries a desolvation penalty that "
                "is overestimated when taken separately. It is not only the desolvation."),
        "concap": ("Concap correction: {n} peptide bonds between pocket residues; the peptide group shared by two neighbouring "
                   "fragments is subtracted once so it is not counted twice. Sum of residues uncorrected -> corrected: {txt}"),
        "fallback": ("WARNING: {n} xtb calculation(s) did not converge with the default settings and were repeated with a higher "
                     "electronic temperature (typically tenths of kcal/mol): {txt}"),
        "concap_off": "No concap correction (--no-concap): the peptide groups shared by neighbours are counted in two residues.",
        "summary": "SUMMARY OF {n} POSES", "rank": "sorted by E_bind of {m}", "pose_col": "pose",
        "res_mean": "Mean contribution per residue ({m}, kcal/mol)", "contact": "contact", "mean": "mean",
        "median": "median", "sd": "sd", "best": "best dist.",
        "filtered": "(the {k} poses with E_bind <= {x:.1f}; the others are left out of the per-residue average)",
        "files": "Files written to", "failed": "Poses that failed",
        "noncommon": "Note: different pockets per pose; energies of different poses are not strictly comparable.",
    },
}


def _f(v, nd=3):
    return "" if v is None else f"{v:.{nd}f}"


def _terms(prefix, d):
    return {f"{prefix}_{t}": _f(None if d is None else d.get(t)) for t in TERMS}


def pose_row(name, r, methods):
    row = {"pose": name, "pocket_residues": r["pocket"]["residues"], "pocket_atoms": r["pocket"]["atoms"],
           "pocket_charge": r["pocket"]["charge"], "ligand_charge": r["ligand"]["charge"],
           "solvent": r["solvent"], "seconds": f"{r['seconds']:.1f}",
           "scf_fallbacks": "; ".join(f"{k} etemp {t}" for k, t in r.get("scf_fallbacks", []))}
    for m in methods:
        d = r["methods"][m]
        row.update({f"{m}_e_bind": _f(d["e_bind"]), f"{m}_e_gas": _f(d["e_gas"]), f"{m}_d_solv": _f(d["d_solv"]),
                    f"{m}_pair_sum": _f(d.get("pair_sum")), f"{m}_residual": _f(d.get("residual")),
                    f"{m}_pair_sum_raw": _f(d.get("pair_sum_raw")), f"{m}_concap_sum": _f(d.get("concap_sum")),
                    f"{m}_pair_sum_gas": _f(d.get("pair_sum_gas")), f"{m}_residual_gas": _f(d.get("residual_gas")),
                    f"{m}_pair_dsolv_sum": _f(d.get("pair_dsolv_sum")), f"{m}_residual_solv": _f(d.get("residual_solv"))})
        row.update(_terms(f"{m}_solv", d["comp"]))
        row.update(_terms(f"{m}_gas", d["comp_gas"]))
    return row


def residue_rows(name, r, methods):
    rows = []
    for w in r["residues"]:
        row = {"pose": name, "residue": w["residue"], "resname": w["resname"], "resnum": w["resnum"], "index": w["index"],
               "in_contact": int(w["contact"]), "min_dist": _f(w["min_dist"]), "ligand_atom": w["lig_atom_name"],
               "ligand_atom_index": w["lig_atom"], "receptor_atom": w["rec_atom"], "min_dist_heavy": _f(w["min_dist_heavy"]),
               "heavy_ligand_atom": w["lig_heavy_name"], "heavy_receptor_atom": w["rec_heavy_atom"],
               "n_atoms_within_cutoff": w["n_atoms_close"]}
        for m in methods:
            d = w[m]
            row.update({f"{m}_pair": _f(d["pair"]), f"{m}_pair_raw": _f(d["pair_raw"]), f"{m}_pair_gas": _f(d["pair_gas"]),
                        f"{m}_dsolv": _f(d["dsolv"])})
            row.update(_terms(f"{m}_solv", d["comp"]))
            row.update(_terms(f"{m}_gas", d["comp_gas"]))
        rows.append(row)
    return rows


def summarize_residues(results, methods, max_ebind=None):
    """Residue statistics over poses for each method: contacts, mean / population sd / median of the pair energy.
    A residue absent from a pose (per-pose pockets) counts as 0. ``max_ebind`` drops the
    poses whose E_bind (first method) is above it."""
    if max_ebind is not None:
        results = [(n, r) for n, r in results if r["methods"][methods[0]]["e_bind"] <= max_ebind]
    n = len(results)
    if n == 0:
        return []
    acc = {}
    for _, r in results:
        for w in r["residues"]:
            a = acc.setdefault(w["residue"], {"resname": w["resname"], "resnum": w["resnum"], "d": [], "c": 0,
                                              **{m: [] for m in methods}})
            a["d"].append(w["min_dist"])
            a["c"] += 1 if w["contact"] else 0
            for m in methods:
                a[m].append(w[m]["pair"])
    out = []

    def stats(vals):
        vals = vals + [0.0] * (n - len(vals))
        mean = sum(vals) / n
        return mean, (sum((v - mean) ** 2 for v in vals) / n) ** 0.5, statistics.median(vals)

    for res, a in acc.items():
        row = {"residue": res, "resname": a["resname"], "resnum": a["resnum"], "n_poses": n, "n_contact": a["c"],
               "contact_fraction": a["c"] / n, "min_dist_mean": sum(a["d"]) / len(a["d"]), "min_dist_best": min(a["d"])}
        for m in methods:
            row[f"{m}_mean"], row[f"{m}_sd"], row[f"{m}_median"] = stats(a[m])
        out.append(row)
    return sorted(out, key=lambda x: x[f"{methods[0]}_mean"])


def _write_csv(path, rows, columns=None):
    if not rows:
        return
    cols = columns or list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(rows)


def write_outputs(outdir, results, methods, max_ebind=None, log_text=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _write_csv(outdir / "poses.csv", [pose_row(n, r, methods) for n, r in results])
    res = [row for n, r in results for row in residue_rows(n, r, methods)]
    if res:
        _write_csv(outdir / "residues.csv", res)
    if len(results) > 1 and res:
        summ = summarize_residues(results, methods, max_ebind)
        if summ:
            _write_csv(outdir / "residues_summary.csv",
                       [{k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in s.items()} for s in summ])
    if log_text is not None:
        (outdir / "qmeku.log").write_text(log_text)


def _table(rows, header, widths, fmt_cols=None):
    lines = ["  ".join(h.rjust(w) if i else h.ljust(w) for i, (h, w) in enumerate(zip(header, widths)))]
    lines.append("  ".join("-" * w for w in widths))
    for r in rows:
        lines.append("  ".join(str(c).rjust(w) if i else str(c).ljust(w) for i, (c, w) in enumerate(zip(r, widths))))
    return lines


def _n(v, w=8, nd=2):
    return " " * (w - 1) + "-" if v is None else f"{v:{w}.{nd}f}"


def format_pose(name, r, methods, lang="en", top=20):
    T = TEXT[lang]
    L = []
    bar = "-" * 100
    L += [bar, f" {T['pose']}: {name}   ({T['ligand_charge']} {r['ligand']['charge']:+d}; "
          + T["pocket_line"].format(r=r["pocket"]["residues"], a=r["pocket"]["atoms"], q=r["pocket"]["charge"]) + f"; {r['seconds']:.0f} s)", bar]
    if r.get("scf_fallbacks"):
        L += [" " + T["fallback"].format(n=len(r["scf_fallbacks"]), txt="; ".join(f"{k} ({t} K)" for k, t in r["scf_fallbacks"])), ""]
    solv_on = r.get("solvent_on", True)
    L += [" " + (T["bind"] if solv_on else T["bind_gas"])]
    for m in methods:
        d = r["methods"][m]
        solv = ("ALPB" if m == "xtb" else "COSMO") if solv_on else "gas"
        L += [f"     {LABEL[m]} ({solv}):".ljust(34) + f"{d['e_bind']:9.2f}"]
    if r["gas"]:
        L += ["", " " + T["balance"]]
        hdr = [""] + [LABEL[m] for m in methods]
        rows = [[T["gas"]] + [f"{r['methods'][m]['e_gas']:.2f}" for m in methods],
                [T["dsolv"]] + [f"{r['methods'][m]['d_solv']:+.2f}" for m in methods],
                [T["ebind"]] + [f"{r['methods'][m]['e_bind']:.2f}" for m in methods]]
        L += ["   " + s for s in _table(rows, hdr, [58] + [11] * len(methods))]
    if "xtb" in methods:
        cg = r["methods"]["xtb"].get("comp_gas")
        cs = r["methods"]["xtb"].get("comp")
        if cg:
            L += ["", f" {T['terms_gas']} - GFN2-xTB"]
            for key, lab in (("es", "t_es"), ("disp", "t_disp"), ("rep", "t_rep"), ("xc", "t_xc"), ("rest", "t_rest")):
                L += [f"     {T[lab]}:".ljust(34) + f"{cg[key]:9.2f}"]
            L += [f"     {T['t_total']}:".ljust(34) + f"{r['methods']['xtb']['e_gas']:9.2f}"]
            if not solv_on:
                L += ["   " + T["terms_note"]]
        if solv_on and cs:
            L += ["", f" {T['terms']} - GFN2-xTB"]
            for key, lab in (("es", "t_es"), ("disp", "t_disp"), ("rep", "t_rep"), ("xc", "t_xc"),
                             ("gsolv", "t_gsolv"), ("rest", "t_rest")):
                L += [f"     {T[lab]}:".ljust(34) + f"{cs[key]:9.2f}"]
            L += [f"     {T['t_total']}:".ljust(34) + f"{r['methods']['xtb']['e_bind']:9.2f}"]
            L += ["   " + T["terms_note"]]
    rows = r["residues"]
    if rows:
        L += _residue_tables(r, methods, lang, top)
    return L


def _residue_tables(r, methods, lang, top):
    """Summary table (energy per method and distances) and, with the gas phase calculated, a second table with the
    gas / solvation split and the xtb terms. Both finish with the sum of the residues, the non-additive remainder and
    the total, so that they reconcile with E_bind."""
    T = TEXT[lang]
    rows = r["residues"]
    main = methods[0]
    order = sorted(rows, key=lambda w: w[main]["pair"])
    shown = order[:top] if top else order
    out = []

    def render(header, widths, body, foot):
        lines = ["  ".join(h.rjust(w) if i else h.ljust(w) for i, (h, w) in enumerate(zip(header, widths))),
                 "  ".join("-" * w for w in widths)]
        lines += ["  ".join(str(c).rjust(w) if i else str(c).ljust(w) for i, (c, w) in enumerate(zip(line, widths)))
                  for line in body]
        lines += ["  ".join("=" * w for w in widths)]
        lines += ["  ".join(str(c).rjust(w) if i else str(c).ljust(w) for i, (c, w) in enumerate(zip(line, widths)))
                  for line in foot]
        return ["   " + x for x in lines]

    # ---- table 1: summary
    out += ["", " " + T["perres"]]
    header = [T["res"], "d (Å)", T["atoms"], "d " + T["heavy"], T["atoms"] + " " + T["heavy"]] + [LABEL[m] for m in methods]
    widths = [18, 6, 15, 8, 15] + [11] * len(methods)
    body = []
    for w in shown:
        line = [w["residue"], f"{w['min_dist']:.2f}", f"{w['lig_atom_name']}-{w['rec_atom']}", f"{w['min_dist_heavy']:.2f}",
                f"{w['lig_heavy_name']}-{w['rec_heavy_atom']}"] + [f"{w[m]['pair']:.2f}" for m in methods]
        body.append(line)
    blank = ["", "", "", ""]
    foot = []
    for lab, key in ((T["sum"], "pair_sum"), (T["rest"], "residual"), ("E_bind", "e_bind")):
        line = [lab] + blank + [f"{r['methods'][m][key]:.2f}" for m in methods]
        foot.append(line)
    out += render(header, widths, body, foot)
    if len(order) > len(shown):
        out += ["   " + T["more"].format(n=len(order) - len(shown))]

    # ---- table 2: gas phase / solvation split and xtb terms
    if r["gas"]:
        out += ["", " " + T["split"]]
        header = [T["res"]]
        widths = [18]
        for m in methods:
            header += [LABEL[m] + " gas", "ΔSolv"]
            widths += [13, 9]
        if "xtb" in methods:
            header += ["ES (gas)", "disp (gas)", "rep (gas)"]
            widths += [9, 10, 9]
        body = []
        for w in shown:
            line = [w["residue"]]
            for m in methods:
                line += [f"{w[m]['pair_gas']:.2f}", f"{w[m]['dsolv']:+.2f}"]
            if "xtb" in methods:
                cg = w["xtb"]["comp_gas"]
                line += [f"{cg['es']:.2f}", f"{cg['disp']:.2f}", f"{cg['rep']:.2f}"]
            body.append(line)
        foot = []
        for lab, gk, sk, term in ((T["sum"], "pair_sum_gas", "pair_dsolv_sum", None), (T["rest"], "residual_gas", "residual_solv", "rest"),
                                  ("E_int gas / ΔSolv", "e_gas", "d_solv", "tot")):
            line = [lab]
            for m in methods:
                d = r["methods"][m]
                line += [f"{d[gk]:.2f}", f"{d[sk]:+.2f}"]
            if "xtb" in methods:
                cg = r["methods"]["xtb"]["comp_gas"]
                sums = {k: sum(w["xtb"]["comp_gas"][k] for w in rows) for k in ("es", "disp", "rep")}
                if term is None:
                    line += [f"{sums[k]:.2f}" for k in ("es", "disp", "rep")]
                elif term == "rest":
                    line += [f"{cg[k] - sums[k]:.2f}" for k in ("es", "disp", "rep")]
                else:
                    line += [f"{cg[k]:.2f}" for k in ("es", "disp", "rep")]
            foot.append(line)
        out += render(header, widths, body, foot)
    out += ["", "   " + T["why"]]
    if r.get("concap"):
        txt = "; ".join(f"{LABEL[m]} {r['methods'][m]['pair_sum_raw']:.2f} → {r['methods'][m]['pair_sum']:.2f}" for m in methods)
        out += ["   " + T["concap"].format(n=r["n_concaps"], txt=txt)]
    else:
        out += ["   " + T["concap_off"]]
    return out


def format_log(results, methods, meta, lang="en", top=20, max_ebind=None, failed=()):
    T = TEXT[lang]
    L = ["=" * 100, " " + T["title"].format(ver=meta["version"]), "=" * 100,
         f" {T['date']}: {datetime.datetime.now():%Y-%m-%d %H:%M}",
         f" {T['receptor']}: {meta['receptor']}", f" {T['poses']}: {meta['poses']} ({len(results)})"]
    ms = []
    for m in methods:
        v = meta["versions"].get(m)
        ms.append(f"{LABEL[m]} ({'xtb' if m == 'xtb' else 'MOPAC'} {v})" + (f", {'ALPB' if m == 'xtb' else 'COSMO'} {meta['solvent']}" if meta["solvent_on"] else ", gas"))
    L += [f" {T['methods']}: " + "; ".join(ms)]
    pocket = (T["pocket_common"] if meta["pocket"] == "common" else T["pocket_each"]).format(c=meta["cutoff"])
    L += [f" {T['pocket']}: {pocket}"]
    L += ["", " " + T["warn"], ""]
    for name, r in results:
        L += format_pose(name, r, methods, lang, top) + [""]
    if len(results) > 1:
        main = methods[0]
        L += ["=" * 100, " " + T["summary"].format(n=len(results)) + "  " + T["rank"].format(m=LABEL[main]), "=" * 100]
        hdr = [T["pose_col"]] + [f"E_bind {LABEL[m]}" for m in methods] + (["E_gas " + LABEL[main], "ΔSolv"] if results[0][1]["gas"] else [])
        widths = [16] + [16] * len(methods) + ([14, 9] if results[0][1]["gas"] else [])
        rows = []
        for n, r in sorted(results, key=lambda x: x[1]["methods"][main]["e_bind"]):
            line = [n] + [f"{r['methods'][m]['e_bind']:.2f}" for m in methods]
            if r["gas"]:
                line += [f"{r['methods'][main]['e_gas']:.2f}", f"{r['methods'][main]['d_solv']:+.2f}"]
            rows.append(line)
        L += ["   " + s for s in _table(rows, hdr, widths)]
        summ = summarize_residues(results, methods, max_ebind)
        if summ and results[0][1]["residues"]:
            k = summ[0]["n_poses"]
            L += ["", " " + T["res_mean"].format(m=LABEL[main]) + ("  " + T["filtered"].format(k=k, x=max_ebind) if k != len(results) else "")]
            hdr = [T["res"], T["contact"], T["mean"], T["median"], T["sd"], T["best"]]
            widths = [12, 9, 8, 8, 7, 11]
            rows = [[s["residue"], f"{s['n_contact']}/{s['n_poses']}", f"{s[main + '_mean']:.2f}", f"{s[main + '_median']:.2f}",
                     f"{s[main + '_sd']:.2f}", f"{s['min_dist_best']:.2f}"] for s in (summ[:top] if top else summ)]
            L += ["   " + s for s in _table(rows, hdr, widths)]
        if meta["pocket"] != "common":
            L += ["", " " + T["noncommon"]]
    for n, e in failed:
        L += ["", f" {T['failed']}: {n}: {e}"]
    return "\n".join(L) + "\n"
