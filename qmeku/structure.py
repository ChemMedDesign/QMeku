"""Pocket models for the quantum-chemical calculations of QMeku.

From a complex PDB (receptor with hydrogens and AMBER names, ligand as HETATM residue) this builds
small systems that xtb / MOPAC can handle:

* :func:`select_residues`   residues with any atom within ``cutoff`` of the ligand;
* :func:`pocket_atoms`      the union of the selected residues (plus single-residue gaps between them), with
                            the cut peptide bonds capped;
* :func:`fragment_atoms`    a single residue capped on both sides (for ligand-residue pair energies).

Capping scheme (hydrogen link atoms, no geometry optimisation):

* the previous residue contributes its C and O, and its C-alpha is replaced by a hydrogen on the C-CA
  bond at 1.09 A (a formyl-like cap);
* the next residue contributes its N (and its H), and its C-alpha (and C-delta for proline) is replaced
  by a hydrogen on the N-C bond at 1.01 A (an amide-like cap);
* a disulfide bridge to a residue that is not included is capped with a hydrogen on the S (1.34 A).

All atoms keep their coordinates of the pose, so distances to the ligand are unchanged. The net charge
is counted from the residue names (ASP/GLU/CYM -1, LYS/ARG/HIP +1, N-terminal NH3+ +1, C-terminal OXT -1);
caps are neutral. Prototype limits: protein residues only, no alternate locations, one chain break rule
(peptide bond when C-N < 1.8 A), no check of the electron parity beyond :func:`n_electrons`.
"""
import math
from collections import OrderedDict
from pathlib import Path

Z = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "NA": 11, "MG": 12, "SI": 14, "P": 15, "S": 16, "CL": 17,
     "K": 19, "CA": 20, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29, "ZN": 30, "SE": 34, "BR": 35, "I": 53}
LINK = {"C": 1.09, "N": 1.01, "S": 1.34}
RES_CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1, "HIP": 1, "CYM": -1}


def _atomic_number(el):
    sym = el.strip().upper()
    if sym in Z:
        return Z[sym]
    try:
        from rdkit import Chem
        num = Chem.GetPeriodicTable().GetAtomicNumber(sym.capitalize())
        if num > 0:
            return num
    except Exception:
        pass
    raise ValueError(f"unknown element {el!r}")


def element(name, column=""):
    """Element symbol of an atom: the PDB element column (77-78) when it is there, else the first letter of the
    atom name (right for protein atoms; wrong for two-letter elements such as Cl or Br, which is why ligands are
    written with the element column). Returned capitalised (``Cl``, ``Br``, ``C``)."""
    col = column.strip()
    if col.isalpha():
        return col.capitalize()
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    raise ValueError("no element in " + repr(name))


class Complex:
    """Atoms of a complex PDB grouped by residue (file order). ``ligand`` is the last residue, or the
    HETATM residue named ``ligand_resname``."""

    def __init__(self, pdb_path, ligand_resname=None):
        self.residues = OrderedDict()           # key -> dict(name, resnum, chain, atoms{name: (el, xyz)})
        order = []
        for line in Path(pdb_path).read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                if line.startswith("ENDMDL"):
                    break
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resnum = line[22:27].strip()
            key = (chain, resnum, resname, line.startswith("HETATM"))
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            digits = resnum.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            r = self.residues.setdefault(key, {"resname": resname, "resnum": int(digits), "resid": resnum,
                                               "chain": chain, "hetatm": key[3], "atoms": OrderedDict()})
            r["atoms"][name] = (element(name, line[76:78]), xyz)
            if key not in order:
                order.append(key)
        self.keys = order
        lig = [k for k in order if k[3] and (k[2] == ligand_resname if ligand_resname else True)]
        if len(lig) != 1:
            raise ValueError(f"expected one ligand residue, found {len(lig)}")
        self.ligand_key = lig[0]
        self.protein_keys = [k for k in order if k != self.ligand_key]
        self._prot_idx = {k: i for i, k in enumerate(self.protein_keys)}
        self.multi_chain = len({self.residues[k]["chain"] for k in self.protein_keys}) > 1

        # Pre-calculate disulfide bridges (S-S) symmetrically for thiol-free cysteines
        self.disulfides = {}
        sg_items = [(k, self.residues[k]["atoms"]["SG"][1]) for k in self.protein_keys
                    if self.residues[k]["resname"] in ("CYX", "CYS")
                    and "SG" in self.residues[k]["atoms"] and "HG" not in self.residues[k]["atoms"]]
        candidates = []
        for i in range(len(sg_items)):
            k1, sg1 = sg_items[i]
            for j in range(i + 1, len(sg_items)):
                k2, sg2 = sg_items[j]
                d = math.dist(sg1, sg2)
                if d < 2.5:
                    candidates.append((d, k1, k2))
        for _, k1, k2 in sorted(candidates):
            if k1 not in self.disulfides and k2 not in self.disulfides:
                self.disulfides[k1] = k2
                self.disulfides[k2] = k1

    def index(self, key):
        """Residue index: position in the file, 1-based."""
        return self.keys.index(key) + 1

    def label(self, key):
        r = self.residues[key]
        base = f"{r['resname']}{r['resid']}"
        if self.multi_chain:
            ch = r["chain"].strip() or "_"
            return f"{ch}:{base}"
        return base

    def ligand_atoms(self):
        r = self.residues[self.ligand_key]
        return [(el, xyz, "lig:" + n) for n, (el, xyz) in r["atoms"].items()]


def _dist(a, b):
    return math.dist(a, b)


def _link(from_xyz, to_xyz, length):
    v = [t - f for f, t in zip(from_xyz, to_xyz)]
    n = math.sqrt(sum(c * c for c in v))
    return tuple(f + c / n * length for f, c in zip(from_xyz, v))


def select_residues(cx, cutoff=4.0):
    """Protein residues with any atom within ``cutoff`` A of any ligand atom."""
    lig = [xyz for _, xyz, _ in cx.ligand_atoms()]
    out = []
    for k in cx.protein_keys:
        close = any(_dist(xyz, l) <= cutoff for _, (el, xyz) in cx.residues[k]["atoms"].items() for l in lig)
        if close:
            out.append(k)
    return out


def _neighbour(cx, key, step):
    """Previous (step=-1) or next (step=+1) residue if it is bonded by a peptide bond, else None."""
    i = cx._prot_idx[key]
    j = i + step
    if j < 0 or j >= len(cx.protein_keys):
        return None
    other = cx.protein_keys[j]
    if cx.residues[key]["chain"] != cx.residues[other]["chain"]:
        return None
    a, b = (other, key) if step < 0 else (key, other)
    ca = cx.residues[a]["atoms"].get("C")
    nb = cx.residues[b]["atoms"].get("N")
    if ca and nb and _dist(ca[1], nb[1]) < 1.8:
        return other
    return None


def _residue_atoms(cx, key, include_set):
    """Atoms of residue ``key`` plus caps for the cut bonds. ``include_set`` are the residues present in
    the system (their bonds are not cut)."""
    res = cx.residues[key]
    atoms = [(el, xyz, f"{cx.label(key)}:{n}") for n, (el, xyz) in res["atoms"].items()]
    # previous residue cut
    prev = _neighbour(cx, key, -1)
    if prev is not None and prev not in include_set:
        p = cx.residues[prev]["atoms"]
        for n in ("C", "O"):
            if n in p:
                atoms.append((p[n][0], p[n][1], f"cap:{cx.label(prev)}:{n}"))
        if "CA" in p and "C" in p:
            atoms.append(("H", _link(p["C"][1], p["CA"][1], LINK["C"]), f"link:{cx.label(prev)}:CA"))
    nxt = _neighbour(cx, key, +1)
    if nxt is not None and nxt not in include_set:
        q = cx.residues[nxt]["atoms"]
        for n in ("N", "H"):
            if n in q:
                atoms.append((q[n][0], q[n][1], f"cap:{cx.label(nxt)}:{n}"))
        for n in (("CA", "CD") if cx.residues[nxt]["resname"] == "PRO" else ("CA",)):
            if n in q and "N" in q:
                atoms.append(("H", _link(q["N"][1], q[n][1], LINK["N"]), f"link:{cx.label(nxt)}:{n}"))
    # disulfide bridge cut (pre-calculated in Complex.__init__)
    if "SG" in res["atoms"] and key in cx.disulfides:
        other = cx.disulfides[key]
        if other not in include_set:
            sg = res["atoms"]["SG"][1]
            sg2 = cx.residues[other]["atoms"]["SG"][1]
            atoms.append(("H", _link(sg, sg2, LINK["S"]), f"link:{cx.label(other)}:SG"))
    return atoms


def bridge_gaps(cx, selected):
    """Add every single residue that lies between two selected residues (peptide-bonded to both): if it
    were left out, both neighbours would put a link atom on the same C-alpha and the two would clash."""
    sel = set(selected)
    for k in cx.protein_keys:
        if k in sel:
            continue
        prev, nxt = _neighbour(cx, k, -1), _neighbour(cx, k, +1)
        if prev in sel and nxt in sel and prev is not None and nxt is not None:
            sel.add(k)
    return [k for k in cx.protein_keys if k in sel]


def pocket_atoms(cx, selected, bridge=True):
    """Atoms (element, xyz, label) of the selected residues with capped cuts, and the residues actually
    included (the selection plus the bridging residues of :func:`bridge_gaps`; ``bridge=False`` keeps
    exactly ``selected``, used to build a pocket without one residue)."""
    used = bridge_gaps(cx, selected) if bridge else [k for k in cx.protein_keys if k in set(selected)]
    inc = set(used)
    out = []
    for k in used:
        out += _residue_atoms(cx, k, inc)
    return out, used


def fragment_atoms(cx, key):
    """A single residue capped on both sides."""
    return _residue_atoms(cx, key, {key})


def peptide_bonds(cx, keys):
    """Pairs (k, k2) of consecutive residues of ``keys`` joined by a peptide bond (C(k)-N(k2) < 1.8 A)."""
    kset = set(keys)
    out = []
    for k in cx.protein_keys:
        if k in kset:
            k2 = _neighbour(cx, k, +1)
            if k2 is not None and k2 in kset:
                out.append((k, k2))
    return out


def concap_atoms(cx, k, k2):
    """The peptide group shared by two consecutive residues, as a small neutral molecule (formamide-like): C and O
    of ``k``, N and H of ``k2``, and a hydrogen link atom in place of each C-alpha (C-delta of proline as well).

    Each of the two capped residue fragments contains this group (as the cap of the other), so a quantity summed
    over the fragments counts it twice; subtracting it once per peptide bond is the standard correction of the
    molecular fragmentation with conjugate caps (MFCC). Returns None if the residues miss the needed atoms."""
    a, b = cx.residues[k]["atoms"], cx.residues[k2]["atoms"]
    if not all(n in a for n in ("C", "O", "CA")) or not all(n in b for n in ("N", "CA")):
        return None
    atoms = [(a["C"][0], a["C"][1], "concap:C"), (a["O"][0], a["O"][1], "concap:O"),
             ("H", _link(a["C"][1], a["CA"][1], LINK["C"]), "concap:linkC_CA"),
             (b["N"][0], b["N"][1], "concap:N")]
    if "H" in b:
        atoms.append((b["H"][0], b["H"][1], "concap:H"))
    for n in (("CA", "CD") if cx.residues[k2]["resname"] == "PRO" else ("CA",)):
        if n in b:
            atoms.append(("H", _link(b["N"][1], b[n][1], LINK["N"]), f"concap:linkN_{n}"))
    return atoms


def residue_charge(cx, key):
    res = cx.residues[key]
    rname = res["resname"]
    atoms = res["atoms"]
    q = RES_CHARGE.get(rname, 0)
    if rname != "ACE" and _neighbour(cx, key, -1) is None and "N" in atoms:
        n_xyz = atoms["N"][1]
        nh = sum(1 for el, xyz in atoms.values() if el == "H" and _dist(n_xyz, xyz) < 1.35)
        req = 2 if rname == "PRO" else 3
        if nh >= req or (nh == 0 and "H3" in atoms):
            q += 1
    if rname != "NME" and _neighbour(cx, key, +1) is None and "OXT" in atoms:
        oxt_xyz = atoms["OXT"][1]
        o_xyz = atoms.get("O", (None, oxt_xyz))[1]
        protonated_cooh = any(el == "H" and min(_dist(oxt_xyz, xyz), _dist(o_xyz, xyz)) < 1.35
                              for el, xyz in atoms.values())
        if not protonated_cooh:
            q -= 1
    return q


def system_charge(cx, keys, ligand_charge=0):
    return sum(residue_charge(cx, k) for k in keys) + ligand_charge


def n_electrons(atoms, charge):
    return sum(_atomic_number(el) for el, _, _ in atoms) - charge


def link_atoms_distance_check(atoms):
    """Smallest distance between any two atoms of the system (a sanity check for overlapping links)."""
    m = 9.9
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            m = min(m, _dist(atoms[i][1], atoms[j][1]))
    return m
