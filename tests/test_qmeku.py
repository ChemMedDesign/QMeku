"""Tests of QMeku. Run from the project root:  python3 -m unittest discover -s tests -v

Needs RDKit. The tests that run xtb or MOPAC are skipped when the program is missing.
"""
import math
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

from qmeku import analysis, cli, engines, ligand, report, structure  # noqa: E402

DATA = ROOT / "examples" / "3ptb"
RECEPTOR = DATA / "receptor.pdb"
HAVE = engines.available()


def crystal_mol():
    tmpl = Chem.MolFromSmiles("NC(=[NH2+])c1ccccc1")
    raw = Chem.MolFromPDBFile(str(DATA / "ligand_crystal.pdb"), removeHs=False, proximityBonding=True)
    return ligand.prepare_ligand(AllChem.AssignBondOrdersFromTemplate(tmpl, Chem.RemoveHs(raw)))


def crystal_complex(tmp):
    text = ligand.build_complex(str(RECEPTOR), crystal_mol())
    path = Path(tmp) / "complex.pdb"
    path.write_text(text)
    return path, text


class TestLigandAndComplex(unittest.TestCase):
    def test_complex(self):
        text = ligand.build_complex(str(RECEPTOR), crystal_mol())
        self.assertEqual(text.count("ENDMDL"), 1)
        lines = [l for l in text.splitlines() if l.startswith("HETATM")]
        self.assertEqual(len(lines), 18)
        self.assertEqual(len({l[12:16] for l in lines}), 18)          # unique atom names

    def test_read_poses_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = crystal_mol()
            m.SetProp("_Name", "mypose")
            w = Chem.SDWriter(str(Path(tmp) / "p.sdf"))
            w.write(m)
            w.write(crystal_mol())
            w.close()
            poses = ligand.read_poses(Path(tmp) / "p.sdf")
        self.assertEqual([n for n, _ in poses], ["mypose", "pose2"])


class TestElements(unittest.TestCase):
    def test_two_letter_elements_are_read_from_the_element_column(self):
        # atom names CL1 / BR1 would be read as C / B if the element came from the name
        m = ligand.prepare_ligand(Chem.MolFromSmiles("Clc1ccc(Br)cc1"), embed=True)
        text = ligand.build_complex(str(RECEPTOR), m)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c.pdb"
            path.write_text(text)
            lig = structure.Complex(str(path), "LIG").ligand_atoms()
        self.assertEqual(sorted({a[0] for a in lig}), ["Br", "C", "Cl", "H"])
        self.assertEqual(structure.n_electrons(lig, 0) % 2, 0)
        self.assertEqual(structure.element("CL1"), "C")                       # name only: the old, wrong behaviour
        self.assertEqual(structure.element("CL1", "CL"), "Cl")


class TestStructure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        path, _ = crystal_complex(cls.tmp.name)
        cls.cx = structure.Complex(str(path), "LIG")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_pocket(self):
        sel = structure.select_residues(self.cx, 4.0)
        labels = {self.cx.label(k) for k in sel}
        self.assertTrue({"ASP189", "SER190", "CYX191", "GLY219", "CYX220"} <= labels)
        atoms, used = structure.pocket_atoms(self.cx, sel)
        self.assertIn("SER217", {self.cx.label(k) for k in used})        # single-residue gap bridged
        self.assertEqual(structure.system_charge(self.cx, used), -1)
        lig = self.cx.ligand_atoms()
        self.assertEqual(structure.n_electrons(atoms + lig, structure.system_charge(self.cx, used, 1)) % 2, 0)
        self.assertGreater(structure.link_atoms_distance_check(atoms + lig), 0.9)

    def test_concap(self):
        by = {self.cx.label(k): k for k in self.cx.protein_keys}
        pairs = structure.peptide_bonds(self.cx, [by["ASP189"], by["SER190"], by["GLY219"]])
        self.assertEqual([(self.cx.label(a), self.cx.label(b)) for a, b in pairs], [("ASP189", "SER190")])
        atoms = structure.concap_atoms(self.cx, by["ASP189"], by["SER190"])
        names = [lab.split(":", 1)[1] for _, _, lab in atoms]
        self.assertEqual(sorted(n for n in names if not n.startswith("link")), ["C", "H", "N", "O"])
        self.assertEqual(sorted(n for n in names if n.startswith("link")), ["linkC_CA", "linkN_CA"])
        self.assertEqual(structure.n_electrons(atoms, 0) % 2, 0)          # neutral closed-shell formamide-like
        c = dict((lab.split(":", 1)[1], xyz) for _, xyz, lab in atoms)
        self.assertAlmostEqual(math.dist(c["linkC_CA"], c["C"]), 1.09, places=2)     # link H on the carbonyl C
        self.assertAlmostEqual(math.dist(c["linkN_CA"], c["N"]), 1.01, places=2)     # link H on the amide N

    def test_proline_concap_has_two_links_on_nitrogen(self):
        pro = next(k for k in self.cx.protein_keys if self.cx.residues[k]["resname"] == "PRO")
        prev = structure._neighbour(self.cx, pro, -1)
        atoms = structure.concap_atoms(self.cx, prev, pro)
        self.assertEqual(sum(1 for _, _, lab in atoms if lab.startswith("concap:link")), 3)   # on C(prev), and on N in place of CA and CD


@unittest.skipUnless(HAVE["xtb"], "xtb not installed")
class TestXtb(unittest.TestCase):
    W = [("O", (0.0, 0.0, 0.117), "O"), ("H", (0.0, 0.757, -0.469), "H1"), ("H", (0.0, -0.757, -0.469), "H2")]

    def test_components_close_and_gas_has_no_solvation(self):
        e = engines.run_xtb(self.W, 0)
        self.assertAlmostEqual(e.total / engines.EH_TO_KCAL, -5.0848, places=3)
        self.assertAlmostEqual(sum(e.components.values()), e.total, places=6)
        g = engines.run_xtb(self.W, 0, solvent=None)
        self.assertEqual(g.components["gsolv"], 0.0)
        self.assertLess(e.components["gsolv"], 0.0)

    def test_scf_fallback_on_the_asp194_ion_pair(self):
        # isolated Asp194 + benzamidinium in the gas phase: with the coordinates of the crystal pose as read from
        # examples/3ptb/poses20_perturbed.sdf (4 decimals) the SCF does not converge with the defaults. Convergence
        # of this system depends on tiny coordinate differences (the unrounded crystal pose converges at 300 K).
        with tempfile.TemporaryDirectory() as tmp:
            name, mol = ligand.read_poses(DATA / "poses20_perturbed.sdf")[0]
            text = ligand.build_complex(str(RECEPTOR), mol, [0])
            path = Path(tmp) / "complex.pdb"
            path.write_text(text)
            cx = structure.Complex(str(path), "LIG")
            k = next(k for k in cx.protein_keys if cx.label(k) == "ASP194")
            atoms = structure.fragment_atoms(cx, k) + cx.ligand_atoms()
            e = engines.run_xtb(atoms, 0, solvent=None)
        self.assertEqual(e.etemp, 600)
        self.assertTrue(math.isfinite(e.total))


@unittest.skipUnless(HAVE["mopac"], "MOPAC not installed")
class TestMopac(unittest.TestCase):
    def test_water(self):
        w = TestXtb.W
        e = engines.run_mopac(w, 0)
        self.assertAlmostEqual(e.total, -56.099, places=2)
        self.assertIsNone(e.components["es"])
        self.assertLess(e.components["gsolv"], 0.0)
        self.assertEqual(engines.run_mopac(w, 0, eps=None).components["gsolv"], 0.0)


@unittest.skipUnless(HAVE["xtb"] and HAVE["mopac"], "xtb and MOPAC needed")
class TestAnalysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        path, cls.text = crystal_complex(cls.tmp.name)
        cls.res = analysis.analyze_pose(path, 1, methods=("xtb", "mopac"), cutoff=3.0, workers=8)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_energies_in_the_expected_range(self):
        self.assertTrue(-30.0 < self.res["methods"]["xtb"]["e_bind"] < -22.0)
        self.assertTrue(-24.0 < self.res["methods"]["mopac"]["e_bind"] < -12.0)

    def test_balance_identities(self):
        for m, d in self.res["methods"].items():
            self.assertAlmostEqual(d["e_bind"], d["e_gas"] + d["d_solv"], places=9)
            self.assertAlmostEqual(d["e_bind"], d["pair_sum"] + d["residual"], places=9)
            self.assertAlmostEqual(d["residual"], d["residual_gas"] + d["residual_solv"], places=6)
            self.assertAlmostEqual(d["pair_sum"], sum(w[m]["pair"] for w in self.res["residues"]), places=9)
            self.assertAlmostEqual(d["pair_sum"], d["pair_sum_raw"] - d["concap_sum"], places=6)
            for w in self.res["residues"]:
                self.assertAlmostEqual(w[m]["pair"], w[m]["pair_gas"] + w[m]["dsolv"], places=9)

    def test_xtb_components_add_up(self):
        d = self.res["methods"]["xtb"]
        self.assertAlmostEqual(sum(d["comp_gas"].values()), d["e_gas"], places=6)
        self.assertAlmostEqual(sum(d["comp"].values()), d["e_bind"], places=6)
        self.assertLess(d["comp_gas"]["es"], -50.0)                      # ion pair: electrostatics dominates
        self.assertLess(d["comp_gas"]["disp"], 0.0)
        self.assertIsNone(self.res["methods"]["mopac"]["comp_gas"]["es"])

    def test_concap_reduces_the_non_additivity(self):
        d = self.res["methods"]["xtb"]
        self.assertGreater(self.res["n_concaps"], 0)
        self.assertLess(abs(d["residual"]), abs(d["e_bind"] - d["pair_sum_raw"]))
        raw = analysis.analyze_pose(Path(self.tmp.name) / "complex.pdb", 1, methods=("xtb",), cutoff=3.0, workers=8,
                                    concap=False)
        self.assertAlmostEqual(raw["methods"]["xtb"]["pair_sum"], d["pair_sum_raw"], places=6)

    def test_residue_table(self):
        top = min(self.res["residues"], key=lambda w: w["xtb"]["pair"])
        self.assertEqual(top["residue"], "ASP189")
        self.assertLess(top["min_dist"], 2.1)
        self.assertTrue(top["lig_atom_name"].startswith("H"))
        self.assertTrue(top["lig_heavy_name"].startswith("N"))
        self.assertTrue(top["rec_heavy_atom"].startswith("O"))
        self.assertLess(top["min_dist_heavy"], 3.0)
        self.assertEqual(len(self.res["residues"]), self.res["pocket"]["selected"])

    def test_input_checks(self):
        path = Path(self.tmp.name) / "complex.pdb"
        with self.assertRaises(analysis.InputError):
            analysis.analyze_pose(path, 0, cutoff=3.0)                    # odd electrons: benzamidinium is +1
        bad = Path(self.tmp.name) / "his.pdb"
        bad.write_text(path.read_text().replace("ASP A 189", "HIS A 189"))
        with self.assertRaises(analysis.InputError):
            analysis.analyze_pose(bad, 1, cutoff=3.0)


class TestReport(unittest.TestCase):
    @staticmethod
    def fake(name, e_bind, rows):
        return name, {"methods": {"xtb": {"e_bind": e_bind}},
                      "residues": [dict(residue=r, resname=r[:3], resnum=int(r[3:]), min_dist=d, contact=d <= 4.0,
                                        xtb={"pair": e}) for r, e, d in rows]}

    def test_summary_statistics_and_filter(self):
        res = [self.fake("a", -20, [("ASP189", -10, 2.0), ("SER190", -2, 3.0)]),
               self.fake("b", -10, [("ASP189", -6, 2.5), ("SER190", -1, 4.5)]),
               self.fake("c", 500, [("ASP189", 100, 1.0), ("SER190", 50, 1.5)])]
        s = {x["residue"]: x for x in report.summarize_residues(res, ["xtb"], max_ebind=0)}
        self.assertEqual(s["ASP189"]["n_poses"], 2)
        self.assertAlmostEqual(s["ASP189"]["xtb_mean"], -8.0)
        self.assertAlmostEqual(s["ASP189"]["xtb_sd"], 2.0)
        self.assertEqual(s["SER190"]["n_contact"], 1)
        self.assertEqual(report.summarize_residues(res, ["xtb"], max_ebind=-100), [])


@unittest.skipUnless(HAVE["xtb"], "xtb not installed")
class TestCli(unittest.TestCase):
    def run_cli(self, *extra):
        tmp = tempfile.TemporaryDirectory()
        sdf = Path(tmp.name) / "p.sdf"
        m = crystal_mol()
        m.SetProp("_Name", "crystal")
        w = Chem.SDWriter(str(sdf))
        w.write(m)
        w.close()
        out = Path(tmp.name) / "out"
        proc = subprocess.run([sys.executable, str(ROOT / "bin" / "qmeku"), str(RECEPTOR), str(sdf), "--outdir", str(out),
                               "--cutoff", "3.0", "--workers", "8", "--method", "xtb", *extra],
                              capture_output=True, text=True)
        return tmp, out, proc

    def test_outputs_and_log_spanish(self):
        tmp, out, proc = self.run_cli("--lang", "es")
        with tmp:
            self.assertEqual(proc.returncode, 0, proc.stderr[-600:])
            for f in ("poses.csv", "residues.csv", "qmeku.log"):
                self.assertTrue((out / f).exists(), f)
            poses = (out / "poses.csv").read_text().splitlines()
            head, row = poses[0].split(","), poses[1].split(",")
            d = dict(zip(head, row))
            self.assertEqual(d["pose"], "crystal")
            self.assertTrue(-30.0 < float(d["xtb_e_bind"]) < -22.0)
            self.assertAlmostEqual(float(d["xtb_e_bind"]), float(d["xtb_e_gas"]) + float(d["xtb_d_solv"]), places=2)
            for c in ("xtb_gas_es", "xtb_gas_disp", "xtb_solv_gsolv", "xtb_residual", "xtb_concap_sum"):
                self.assertIn(c, head)
            res = (out / "residues.csv").read_text().splitlines()
            rh = res[0].split(",")
            for c in ("xtb_pair", "xtb_pair_raw", "xtb_pair_gas", "xtb_dsolv", "min_dist", "min_dist_heavy", "in_contact"):
                self.assertIn(c, rh)
            log = (out / "qmeku.log").read_text()
            for text in ("ENERGÍA DE UNIÓN APROXIMADA", "Balance de E_bind", "Contribución por residuo",
                         "resto no aditivo", "Corrección de concaps", "ASP189"):
                self.assertIn(text, log)

    def test_default_english_log_and_total_only(self):
        tmp, out, proc = self.run_cli("--total-only", "--no-gas")
        with tmp:
            self.assertEqual(proc.returncode, 0, proc.stderr[-600:])
            log = (out / "qmeku.log").read_text()
            self.assertIn("APPROXIMATE BINDING ENERGY", log)
            self.assertNotIn("Per-residue contribution", log)
            self.assertFalse((out / "residues.csv").exists())


def _pdb_line(rec, serial, name, resname, chain, resseq, x, y, z, el):
    return f"{rec:<6s}{serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resseq:>4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}"


class TestDisulfidesTerminiAndChains(unittest.TestCase):
    def test_disulfide_string_order_symmetry_and_unpaired_cyx(self):
        lines = [
            _pdb_line("ATOM", 1, "N", "CYX", "A", 22, 0.0, 0.0, 0.0, "N"),
            _pdb_line("ATOM", 2, "CA", "CYX", "A", 22, 1.4, 0.0, 0.0, "C"),
            _pdb_line("ATOM", 3, "C", "CYX", "A", 22, 2.5, 0.0, 0.0, "C"),
            _pdb_line("ATOM", 4, "O", "CYX", "A", 22, 2.5, 1.2, 0.0, "O"),
            _pdb_line("ATOM", 5, "SG", "CYX", "A", 22, 0.0, 2.0, 0.0, "S"),
            # CYX 157 has '157' < '22' lexicographically, testing integer index pairing
            _pdb_line("ATOM", 6, "N", "CYX", "A", 157, 0.0, 6.0, 0.0, "N"),
            _pdb_line("ATOM", 7, "CA", "CYX", "A", 157, 1.4, 6.0, 0.0, "C"),
            _pdb_line("ATOM", 8, "C", "CYX", "A", 157, 2.5, 6.0, 0.0, "C"),
            _pdb_line("ATOM", 9, "O", "CYX", "A", 157, 2.5, 7.2, 0.0, "O"),
            _pdb_line("ATOM", 10, "SG", "CYX", "A", 157, 0.0, 4.05, 0.0, "S"),
            # Third CYX slightly further (2.25 A from CYX 157) must NOT overwrite symmetric pair
            _pdb_line("ATOM", 11, "N", "CYX", "A", 158, 0.0, 9.0, 0.0, "N"),
            _pdb_line("ATOM", 12, "CA", "CYX", "A", 158, 1.4, 9.0, 0.0, "C"),
            _pdb_line("ATOM", 13, "C", "CYX", "A", 158, 2.5, 9.0, 0.0, "C"),
            _pdb_line("ATOM", 14, "O", "CYX", "A", 158, 2.5, 10.2, 0.0, "O"),
            _pdb_line("ATOM", 15, "SG", "CYX", "A", 158, 0.0, 6.30, 0.0, "S"),
            # Protonated CYS with HG within 2.1 A of CYX 22 must be ignored as thiol
            _pdb_line("ATOM", 16, "N", "CYS", "A", 200, -4.0, 2.0, 0.0, "N"),
            _pdb_line("ATOM", 17, "CA", "CYS", "A", 200, -3.0, 2.0, 0.0, "C"),
            _pdb_line("ATOM", 18, "C", "CYS", "A", 200, -3.0, 3.5, 0.0, "C"),
            _pdb_line("ATOM", 19, "O", "CYS", "A", 200, -4.0, 4.2, 0.0, "O"),
            _pdb_line("ATOM", 20, "SG", "CYS", "A", 200, -2.1, 2.0, 0.0, "S"),
            _pdb_line("ATOM", 21, "HG", "CYS", "A", 200, -2.1, 3.3, 0.0, "H"),
            _pdb_line("HETATM", 22, "C1", "LIG", "L", 1, 10.0, 10.0, 10.0, "C"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "ss.pdb"
            p.write_text("\n".join(lines) + "\nEND\n")
            cx = structure.Complex(str(p), "LIG")
            by = {cx.label(k): k for k in cx.protein_keys}
            self.assertEqual(cx.disulfides.get(by["CYX22"]), by["CYX157"])
            self.assertEqual(cx.disulfides.get(by["CYX157"]), by["CYX22"])
            self.assertNotIn(by["CYX158"], cx.disulfides)
            self.assertNotIn(by["CYS200"], cx.disulfides)
            with self.assertRaises(analysis.InputError):
                analysis.check_residues(cx, cx.protein_keys)  # unpaired CYX158 rejected

    def test_proline_n_terminus_and_corrupted_internal_termini(self):
        # N-terminal protonated PRO (H2 + H3 -> +1, 20 atoms / closed shell when C-capped)
        lines = [
            _pdb_line("ATOM", 1, "N", "PRO", "A", 1, 0.0, 0.0, 0.0, "N"),
            _pdb_line("ATOM", 2, "H2", "PRO", "A", 1, -0.5, 0.8, 0.0, "H"),
            _pdb_line("ATOM", 3, "H3", "PRO", "A", 1, -0.5, -0.8, 0.0, "H"),
            _pdb_line("ATOM", 4, "CA", "PRO", "A", 1, 1.45, 0.0, 0.0, "C"),
            _pdb_line("ATOM", 5, "HA", "PRO", "A", 1, 1.8, 0.9, 0.5, "H"),
            _pdb_line("ATOM", 6, "CB", "PRO", "A", 1, 1.9, -1.2, 0.7, "C"),
            _pdb_line("ATOM", 7, "HB2", "PRO", "A", 1, 2.9, -1.2, 0.9, "H"),
            _pdb_line("ATOM", 8, "HB3", "PRO", "A", 1, 1.4, -1.3, 1.6, "H"),
            _pdb_line("ATOM", 9, "CG", "PRO", "A", 1, 1.4, -2.3, -0.2, "C"),
            _pdb_line("ATOM", 10, "HG2", "PRO", "A", 1, 2.1, -3.1, -0.2, "H"),
            _pdb_line("ATOM", 11, "HG3", "PRO", "A", 1, 0.4, -2.6, 0.1, "H"),
            _pdb_line("ATOM", 12, "CD", "PRO", "A", 1, 1.2, -1.6, -1.5, "C"),
            _pdb_line("ATOM", 13, "HD2", "PRO", "A", 1, 2.1, -1.5, -2.0, "H"),
            _pdb_line("ATOM", 14, "HD3", "PRO", "A", 1, 0.4, -2.0, -2.1, "H"),
            _pdb_line("ATOM", 15, "C", "PRO", "A", 1, 2.0, 0.1, -1.4, "C"),
            _pdb_line("ATOM", 16, "O", "PRO", "A", 1, 1.3, 0.1, -2.4, "O"),
            # Peptide-bonded GLY 2 (C-N = 1.33 A) with C-terminal OXT (-1)
            _pdb_line("ATOM", 17, "N", "GLY", "A", 2, 3.33, 0.1, -1.4, "N"),
            _pdb_line("ATOM", 18, "H", "GLY", "A", 2, 3.85, 0.1, -0.5, "H"),
            _pdb_line("ATOM", 19, "CA", "GLY", "A", 2, 4.1, 0.1, -2.6, "C"),
            _pdb_line("ATOM", 20, "HA2", "GLY", "A", 2, 3.9, 1.0, -3.2, "H"),
            _pdb_line("ATOM", 21, "HA3", "GLY", "A", 2, 3.9, -0.8, -3.2, "H"),
            _pdb_line("ATOM", 22, "C", "GLY", "A", 2, 5.6, 0.1, -2.3, "C"),
            _pdb_line("ATOM", 23, "O", "GLY", "A", 2, 6.0, 0.1, -1.1, "O"),
            _pdb_line("ATOM", 24, "OXT", "GLY", "A", 2, 6.4, 0.1, -3.3, "O"),
            _pdb_line("HETATM", 25, "C1", "LIG", "L", 1, 2.0, 3.0, 0.0, "C"),
            _pdb_line("HETATM", 26, "H1", "LIG", "L", 1, 2.0, 4.1, 0.0, "H"),
            _pdb_line("HETATM", 27, "H2", "LIG", "L", 1, 3.0, 2.6, 0.0, "H"),
            _pdb_line("HETATM", 28, "H3", "LIG", "L", 1, 1.5, 2.6, 0.9, "H"),
            _pdb_line("HETATM", 29, "H4", "LIG", "L", 1, 1.5, 2.6, -0.9, "H"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "pro.pdb"
            p.write_text("\n".join(lines) + "\nEND\n")
            cx = structure.Complex(str(p), "LIG")
            analysis.check_residues(cx, cx.protein_keys)
            pro_k, gly_k = cx.protein_keys
            self.assertEqual(structure.residue_charge(cx, pro_k), 1)
            self.assertEqual(structure.residue_charge(cx, gly_k), -1)
            self.assertEqual(structure.n_electrons(structure.fragment_atoms(cx, pro_k), 1) % 2, 0)
            self.assertEqual(structure.n_electrons(structure.fragment_atoms(cx, gly_k), -1) % 2, 0)

            # Corrupt PRO 1 with an internal OXT (simulating the 1BCU insertion-code defect)
            bad_lines = lines[:16] + [_pdb_line("ATOM", 99, "OXT", "PRO", "A", 1, 3.2, 0.1, -1.5, "O")] + lines[16:]
            bad_p = Path(tmp) / "bad_oxt.pdb"
            bad_p.write_text("\n".join(bad_lines) + "\nEND\n")
            bad_cx = structure.Complex(str(bad_p), "LIG")
            self.assertEqual(structure.residue_charge(bad_cx, bad_cx.protein_keys[0]), 1)  # internal OXT ignored
            with self.assertRaises(analysis.InputError):
                analysis.check_residues(bad_cx, bad_cx.protein_keys)

    def test_per_residue_parity_under_total_only_and_clash_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, text = crystal_complex(tmp)
            # Remove one H from two different pocket residues so each is odd-electron, while sum is even!
            lines = text.splitlines()
            dropped = set()
            corrupted = []
            for line in lines:
                if line.startswith("ATOM"):
                    res = (line[17:20].strip(), line[22:27].strip())
                    name = line[12:16].strip()
                    if res in (("ASP", "189"), ("SER", "190")) and name == "HA" and res not in dropped:
                        dropped.add(res)
                        continue
                corrupted.append(line)
            odd_pair_pdb = Path(tmp) / "odd_pair.pdb"
            odd_pair_pdb.write_text("\n".join(corrupted) + "\n")
            with self.assertRaisesRegex(analysis.InputError, "odd number of electrons"):
                analysis.analyze_pose(odd_pair_pdb, 1, methods=("xtb",), cutoff=3.0, pairs=False)

            # Place a ligand atom within 0.2 A of ASP189 OD2 to trigger overlap InputError
            clash_lines = []
            od2_xyz = None
            for line in lines:
                if line.startswith("ATOM") and line[17:20].strip() == "ASP" and line[22:27].strip() == "189" and line[12:16].strip() == "OD2":
                    od2_xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                if line.startswith("HETATM") and od2_xyz is not None and line[12:16].strip() == "C1":
                    line = f"{line[:30]}{od2_xyz[0] + 0.10:8.3f}{od2_xyz[1]:8.3f}{od2_xyz[2]:8.3f}{line[54:]}"
                clash_lines.append(line)
            clash_pdb = Path(tmp) / "clash.pdb"
            clash_pdb.write_text("\n".join(clash_lines) + "\n")
            with self.assertRaisesRegex(analysis.InputError, "overlapping atoms"):
                analysis.analyze_pose(clash_pdb, 1, methods=("xtb",), cutoff=3.0, pairs=False)

    def test_multichain_labels_and_homodimer_summary(self):
        lines = [
            _pdb_line("ATOM", 1, "N", "ASP", "A", 189, 0.0, 0.0, 0.0, "N"),
            _pdb_line("ATOM", 2, "CA", "ASP", "A", 189, 1.4, 0.0, 0.0, "C"),
            _pdb_line("ATOM", 3, "C", "ASP", "A", 189, 2.5, 0.0, 0.0, "C"),
            _pdb_line("ATOM", 4, "O", "ASP", "A", 189, 2.5, 1.2, 0.0, "O"),
            _pdb_line("ATOM", 5, "N", "ASP", "B", 189, 0.0, 8.0, 0.0, "N"),
            _pdb_line("ATOM", 6, "CA", "ASP", "B", 189, 1.4, 8.0, 0.0, "C"),
            _pdb_line("ATOM", 7, "C", "ASP", "B", 189, 2.5, 8.0, 0.0, "C"),
            _pdb_line("ATOM", 8, "O", "ASP", "B", 189, 2.5, 9.2, 0.0, "O"),
            _pdb_line("HETATM", 9, "C1", "LIG", "L", 1, 0.0, 4.0, 0.0, "C"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "dimer.pdb"
            p.write_text("\n".join(lines) + "\nEND\n")
            cx = structure.Complex(str(p), "LIG")
            self.assertTrue(cx.multi_chain)
            self.assertEqual([cx.label(k) for k in cx.protein_keys], ["A:ASP189", "B:ASP189"])

        poses = [
            ("p1", {"methods": {"xtb": {"e_bind": -12.0}},
                    "residues": [
                        dict(residue="A:ASP189", resname="ASP", resnum=189, min_dist=2.2, contact=True, xtb={"pair": -8.0}),
                        dict(residue="B:ASP189", resname="ASP", resnum=189, min_dist=2.8, contact=True, xtb={"pair": -4.0}),
                    ]}),
            ("p2", {"methods": {"xtb": {"e_bind": -10.0}},
                    "residues": [
                        dict(residue="A:ASP189", resname="ASP", resnum=189, min_dist=2.4, contact=True, xtb={"pair": -6.0}),
                        dict(residue="B:ASP189", resname="ASP", resnum=189, min_dist=4.5, contact=False, xtb={"pair": -2.0}),
                    ]}),
        ]
        summ = {r["residue"]: r for r in report.summarize_residues(poses, ["xtb"])}
        self.assertEqual(set(summ), {"A:ASP189", "B:ASP189"})
        self.assertAlmostEqual(summ["A:ASP189"]["xtb_mean"], -7.0)
        self.assertAlmostEqual(summ["B:ASP189"]["xtb_mean"], -3.0)
        self.assertLessEqual(summ["A:ASP189"]["contact_fraction"], 1.0)
        self.assertAlmostEqual(summ["B:ASP189"]["contact_fraction"], 0.5)


class TestElementsAndLigandLookup(unittest.TestCase):
    def test_boron_selenium_and_hetatm_ligand_resname(self):
        atoms = [("B", (0.0, 0.0, 0.0), "B1"), ("Se", (2.0, 0.0, 0.0), "SE1"), ("H", (0.0, 1.2, 0.0), "H1")]
        self.assertEqual(structure.n_electrons(atoms, 0), 5 + 34 + 1)
        # Ligand named PRO in HETATM must not collide with protein ATOM PRO residues
        lines = [
            _pdb_line("ATOM", 1, "N", "PRO", "A", 10, 0.0, 0.0, 0.0, "N"),
            _pdb_line("ATOM", 2, "CA", "PRO", "A", 10, 1.4, 0.0, 0.0, "C"),
            _pdb_line("HETATM", 3, "N1", "PRO", "L", 1, 8.0, 8.0, 8.0, "N"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "pro_lig.pdb"
            p.write_text("\n".join(lines) + "\nEND\n")
            cx = structure.Complex(str(p), "PRO")
            self.assertEqual(cx.ligand_key, ("L", "1", "PRO", True))
            self.assertEqual(len(cx.protein_keys), 1)


class TestReportAndGasBookkeeping(unittest.TestCase):
    @staticmethod
    def _mock_energy(atoms, q):
        # Deterministic toy energy where gsolv = 0 in gas phase and components add up to total
        total = -0.25 * len(atoms) - 0.1 * q
        comps = {"es": -0.10 * len(atoms), "xc": -0.02 * len(atoms), "disp": -0.08 * len(atoms),
                 "rep": 0.02 * len(atoms), "gsolv": 0.0}
        comps["rest"] = total - sum(comps.values())
        return engines.Energy(total, comps)

    def test_solvent_none_populates_gas_fields_and_formats_cleanly(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, _ = crystal_complex(tmp)
            with mock.patch("qmeku.analysis.make_runner", return_value=self._mock_energy):
                res = analysis.analyze_pose(path, 1, methods=("xtb",), solvent="none", cutoff=3.0, workers=2)
            d = res["methods"]["xtb"]
            self.assertEqual(d["e_bind"], d["e_gas"])
            self.assertIsNone(d["d_solv"])
            self.assertIsNone(d["comp"])
            self.assertIsNotNone(d["comp_gas"])
            self.assertAlmostEqual(sum(d["comp_gas"].values()), d["e_gas"], places=9)
            self.assertAlmostEqual(d["pair_sum_gas"], d["pair_sum"], places=9)
            self.assertAlmostEqual(d["residual_gas"], d["residual"], places=9)
            for w in res["residues"]:
                self.assertEqual(w["xtb"]["pair_gas"], w["xtb"]["pair"])
                self.assertIsNone(w["xtb"]["dsolv"])
                self.assertIsNone(w["xtb"]["comp"])
                self.assertIsNotNone(w["xtb"]["comp_gas"])
            meta_gas = {"version": "1.0", "receptor": "rec.pdb", "poses": "lig.sdf",
                        "versions": {"xtb": "mock"}, "solvent": "none", "solvent_on": False,
                        "pocket": "per-pose", "cutoff": 3.0}
            log_en = report.format_log([("crystal", res)], ["xtb"], meta_gas, lang="en", top=0)
            self.assertIn("GAS-PHASE INTERACTION ENERGY", log_en)
            self.assertNotIn("Balance of E_bind", log_en)

    def test_no_gas_still_prints_solvent_gfn2_terms_and_top_exceeding_15(self):
        def mock_solv_energy(atoms, q):
            total = -0.30 * len(atoms)
            comps = {"es": -0.10 * len(atoms), "xc": -0.02 * len(atoms), "disp": -0.08 * len(atoms),
                     "rep": 0.02 * len(atoms), "gsolv": -0.05 * len(atoms)}
            comps["rest"] = total - sum(comps.values())
            return engines.Energy(total, comps)

        with tempfile.TemporaryDirectory() as tmp:
            path, _ = crystal_complex(tmp)
            with mock.patch("qmeku.analysis.make_runner", return_value=mock_solv_energy):
                res = analysis.analyze_pose(path, 1, methods=("xtb",), solvent="water", gas=False,
                                            cutoff=4.5, workers=2)
            self.assertGreater(len(res["residues"]), 15)
            poses = [("p1", res), ("p2", res)]
            meta_solv = {"version": "1.0", "receptor": "rec.pdb", "poses": "lig.sdf",
                         "versions": {"xtb": "mock"}, "solvent": "water", "solvent_on": True,
                         "pocket": "common", "cutoff": 4.5}
            log_en = report.format_log(poses, ["xtb"], meta_solv, lang="en", top=20)
            self.assertIn("GFN2-xTB terms in implicit solvent", log_en)
            # Summary table with top=20 must include more than 15 residues
            summary_block = log_en.split("Mean contribution per residue", 1)[1]
            data_lines = [ln for ln in summary_block.splitlines() if ln.strip() and ln.strip()[0].isalpha()
                          and not ln.strip().startswith(("residue", "pose"))]
            self.assertEqual(len(data_lines), min(20, len(res["residues"])))


class TestTimeoutsAndCancellation(unittest.TestCase):
    def test_timeout_validation_and_subprocess_wrapping(self):
        w = [("O", (0.0, 0.0, 0.117), "O"), ("H", (0.0, 0.757, -0.469), "H1"), ("H", (0.0, -0.757, -0.469), "H2")]
        with self.assertRaises(ligand.QmekuError):
            engines.run_xtb(w, 0, timeout=-5)
        with self.assertRaises(ligand.QmekuError):
            engines.run_mopac(w, 0, timeout=0)
        with self.assertRaises(ligand.QmekuError):
            ligand.prepare_receptor("in.pdb", "out.pdb", timeout=-1)

        with mock.patch("qmeku.engines.xtb_bin", return_value="/bin/true"), \
                mock.patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="xtb", timeout=1.5)):
            with self.assertRaisesRegex(ligand.QmekuError, "xtb timed out"):
                engines.run_xtb(w, 0, timeout=1.5)

        with mock.patch("qmeku.engines.mopac_bin", return_value="/bin/true"), \
                mock.patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="mopac", timeout=2.0)):
            with self.assertRaisesRegex(ligand.QmekuError, "mopac timed out"):
                engines.run_mopac(w, 0, timeout=2.0)

        with tempfile.TemporaryDirectory() as tmp:
            in_pdb = Path(tmp) / "in.pdb"
            in_pdb.write_text(_pdb_line("ATOM", 1, "N", "GLY", "A", 1, 0.0, 0.0, 0.0, "N") + "\nEND\n")
            with mock.patch("shutil.which", return_value="/bin/true"), \
                    mock.patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="pdb2pqr", timeout=3.0)):
                with self.assertRaisesRegex(ligand.QmekuError, "pdb2pqr timed out"):
                    ligand.prepare_receptor(in_pdb, Path(tmp) / "out.pdb", timeout=3.0)

    def test_xtb_single_total_deadline_across_scf_retries(self):
        w = [("O", (0.0, 0.0, 0.117), "O"), ("H", (0.0, 0.757, -0.469), "H1"), ("H", (0.0, -0.757, -0.469), "H2")]
        seen_timeouts = []
        clock = [100.0]

        def fake_monotonic():
            return clock[0]

        def fake_run(*args, **kwargs):
            seen_timeouts.append(kwargs.get("timeout"))
            clock[0] += 4.0  # advance 4s out of a 10s total budget on each non-converging attempt
            return subprocess.CompletedProcess(args[0], returncode=1, stdout="SCF did not converge", stderr="")

        with mock.patch("qmeku.engines.xtb_bin", return_value="/bin/true"), \
                mock.patch("qmeku.engines.time.monotonic", side_effect=fake_monotonic), \
                mock.patch("subprocess.run", side_effect=fake_run):
            with self.assertRaisesRegex(ligand.QmekuError, "xtb timed out"):
                engines.run_xtb(w, 0, timeout=10.0)
        # Attempt 1 gets 10s, Attempt 2 gets 6s, Attempt 3 gets 2s, then budget is exhausted before Attempt 4!
        self.assertEqual(seen_timeouts, [10.0, 6.0, 2.0])

    def test_worker_pool_cancels_pending_jobs_on_first_failure(self):
        executed = []

        def failing_runner(atoms, q):
            executed.append(len(atoms))
            if len(executed) == 1:
                raise ligand.QmekuError("synthetic SCF failure")
            time.sleep(0.2)
            return engines.Energy(-1.0, {k: 0.0 for k in engines.COMP_KEYS})

        with tempfile.TemporaryDirectory() as tmp:
            path, _ = crystal_complex(tmp)
            with mock.patch("qmeku.analysis.make_runner", return_value=failing_runner):
                with self.assertRaisesRegex(ligand.QmekuError, "synthetic SCF failure"):
                    analysis.analyze_pose(path, 1, methods=("xtb",), cutoff=4.0, workers=1)
        # With workers=1 and immediate cancellation, queued jobs after the first failure do not all execute
        self.assertLessEqual(len(executed), 2)

    def test_cli_negative_timeout_and_mopac_invalid_solvent(self):
        args = cli.build_parser().parse_args(["rec.pdb", "lig.sdf", "--outdir", "out", "--timeout", "-10"])
        with self.assertRaisesRegex(ligand.QmekuError, "--timeout must be >= 0"):
            cli.run(args)
        with self.assertRaisesRegex(ligand.QmekuError, "MOPAC solvent must be 'water'"):
            engines.make_runner("mopac", True, solvent="methanol")


class TestPymolView(unittest.TestCase):
    def test_color_scale(self):
        from qmeku import pymol_view as pv
        self.assertEqual(pv.energy_color(-10, -10, 2), pv.RED)          # most negative: red
        self.assertEqual(pv.energy_color(0, -10, 2), pv.GREY)           # zero: light grey
        self.assertEqual(pv.energy_color(2, -10, 2), pv.LIGHT_BLUE)     # positive: light blue
        self.assertEqual(pv.energy_color(1, -10, 0), pv.GREY)           # no positive scale: grey

    def test_labels_and_scene(self):
        from qmeku import pymol_view as pv
        self.assertEqual(pv.parse_label("ALA221A"), (None, "ALA", "221A"))
        self.assertEqual(pv.parse_label("B:ASP189"), ("B", "ASP", "189"))
        rows = [{"chain": None, "resname": "ASP", "resi": "189", "e": -9.0, "contact": True, "lig_idx": 9, "rec_atom": "OD2"},
                {"chain": None, "resname": "SER", "resi": "195", "e": 1.0, "contact": True, "lig_idx": None, "rec_atom": None},
                {"chain": None, "resname": "GLY", "resi": "10", "e": -0.5, "contact": False, "lig_idx": None, "rec_atom": None}]
        cmds = pv.scene_commands(rows, "ligand_x.sdf", -9.0, 1.0)
        text = "\n".join(cmds)
        self.assertIn('"Asp189"', text)                                 # capitalised, no energy value
        self.assertNotIn("-9.0", text.split("alter receptor, b=0.0")[0])
        self.assertIn("set sphere_scale, 0.2, ligand", text)
        self.assertIn("set stick_radius, 0.18, ligand", text)
        self.assertIn("color paleyellow, ligand", text)
        self.assertIn("distance qm_polar_1, ligand and index 9", text)
        self.assertNotIn("resi 10)", text.split("select qm_pocket")[1])  # residue out of contact not in the view

    def test_flags(self):
        ap = cli.build_parser()
        a = ap.parse_args(["r.pdb", "p.sdf", "--outdir", "o"])
        self.assertFalse(a.no_pymol)
        self.assertEqual(a.pymol_poses, 3)
        self.assertEqual(a.pymol_method, "both")
        self.assertEqual(ap.parse_args(["r.pdb", "p.sdf", "--outdir", "o", "--pymol-method", "mopac"]).pymol_method, "mopac")
        self.assertTrue(ap.parse_args(["r.pdb", "p.sdf", "--outdir", "o", "--no-pymol"]).no_pymol)


if __name__ == "__main__":
    unittest.main()
