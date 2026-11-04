"""Preparation of ligands (RDKit) and of the complex PDB that the rest of QMeku reads.

The receptor must already have hydrogens and AMBER residue/atom names (``prepare_receptor`` runs PDB2PQR for
that). The ligand enters as an RDKit molecule with bond orders and formal charges, so its net charge is known
without guessing it from distances.
"""
import shutil
import subprocess
import tempfile
from pathlib import Path


class QmekuError(RuntimeError):
    """Problem with the inputs or with an external program."""


def _read_text(x):
    """PDB text, or the content of the file when ``x`` is a path (a string without newlines)."""
    s = str(x)
    return s if "\n" in s else Path(s).read_text()


def prepare_receptor(pdb_in, pdb_out, ph=7.4, timeout=None):
    """Protonate a receptor with PDB2PQR (AMBER names, pH ``ph``) and keep only ATOM records.

    Alternate locations other than ' ' and 'A', HETATM records and waters are dropped."""
    if timeout is not None and timeout <= 0:
        raise QmekuError(f"timeout must be > 0 (got {timeout})")
    exe = shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
    if exe is None:
        raise QmekuError("pdb2pqr not found: pip install pdb2pqr")
    with tempfile.TemporaryDirectory() as tmp:
        heavy = Path(tmp) / "heavy.pdb"
        lines = [l[:16] + " " + l[17:] for l in Path(pdb_in).read_text().splitlines()
                 if l.startswith("ATOM") and len(l) > 16 and l[16] in " A"]
        heavy.write_text("\n".join(lines) + "\nEND\n")
        out = Path(tmp) / "h.pdb"
        try:
            subprocess.run([exe, "--ff=AMBER", "--ffout=AMBER", f"--with-ph={ph}", "--drop-water", "--keep-chain",
                            "--pdb-output", str(out), str(heavy), str(Path(tmp) / "h.pqr")],
                           check=True, timeout=timeout, capture_output=True, text=True)
        except subprocess.TimeoutExpired:
            raise QmekuError(f"pdb2pqr timed out after {timeout} s") from None
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or str(exc)).strip()[-400:]
            raise QmekuError(f"pdb2pqr failed: {detail}") from None
        atoms = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
    Path(pdb_out).write_text("\n".join(atoms) + "\nEND\n")
    return str(pdb_out)


def receptor_atom_lines(receptor_pdb):
    """ATOM records of the first model of a receptor PDB (path or text)."""
    atoms = []
    for line in _read_text(receptor_pdb).splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM"):
            atoms.append(line.rstrip("\n"))
    if not atoms:
        raise QmekuError("no ATOM records in the receptor")
    return atoms


def prepare_ligand(mol, embed=False, seed=7):
    """RDKit molecule with explicit hydrogens and 3D coordinates.

    ``mol`` must carry bond orders and formal charges (SMILES, SDF, MOL2 or a template). Hydrogens are added
    with coordinates; with ``embed`` (or no conformer) a conformer is generated with ETKDG + MMFF, which is
    not a pose in the binding site."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    has_conf = mol.GetNumConformers() > 0
    m = Chem.AddHs(mol, addCoords=has_conf and not embed)
    if embed or not has_conf:
        if AllChem.EmbedMolecule(m, randomSeed=seed) != 0:
            raise QmekuError("could not generate 3D coordinates")
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        except Exception:
            pass
    return m


def read_poses(path):
    """Ligands of an SDF/MOL file as ``[(name, rdkit_mol)]`` (explicit hydrogens added when missing)."""
    from rdkit import Chem
    mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
    if not mols:
        raise QmekuError("no molecules read from " + str(path))
    out = []
    for i, m in enumerate(mols):
        name = m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else f"pose{i + 1}"
        out.append((name, prepare_ligand(m)))
    return out


def ligand_pdb_lines(mol, conf_id=-1, resname="LIG", chain="L", resnum=1):
    """HETATM records of one conformer with unique element-based names (C1, C2, N1, H1 ...)."""
    conf = mol.GetConformer(conf_id)
    counts, lines = {}, []
    for a in mol.GetAtoms():
        el = a.GetSymbol()
        counts[el] = counts.get(el, 0) + 1
        name = f"{el.upper()}{counts[el]}"
        name = (" " + name).ljust(4) if len(name) < 4 else name
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append("HETATM%5d %4s %-3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                     % (0, name, resname[:3], chain, resnum, p.x, p.y, p.z, el.upper()))
    return lines


def build_complex(receptor, ligand, conf_ids=None):
    """Complex as multi-model PDB text (receptor atoms + the ligand as the last residue, one model per
    conformer)."""
    rec = receptor_atom_lines(receptor)
    ids = list(conf_ids) if conf_ids is not None else [c.GetId() for c in ligand.GetConformers()]
    if not ids:
        raise QmekuError("the ligand has no conformers")
    out = []
    for cid in ids:
        serial, model = 0, []
        for line in rec:
            serial += 1
            model.append(line[:6] + "%5d" % serial + line[11:])
        model.append("TER")
        for line in ligand_pdb_lines(ligand, cid):
            serial += 1
            model.append(line[:6] + "%5d" % serial + line[11:])
        out += model + ["ENDMDL"]
    return "\n".join(out) + "\nEND\n"
