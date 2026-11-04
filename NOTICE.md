# Provenance, authorship and licences of QMeku

## Original code

All the code of the `qmeku/` package (pocket selection, closure of the cut bonds, bridging residues, concap
correction, execution of xtb and MOPAC, per-residue analysis, outputs and logging) is original work of the
**SB2CS lab** for QMeku. It is distributed under the **GNU General Public License, version 3 or later** (see `LICENSE`).

## History

QMeku began as an independent tool for analysing docking poses with semiempirical methods. Earlier versions
included, as an optional complement, a classical engine in C (`engine/`) derived from **MM-ISMSA / cMMISMSA**
by Álvaro Cortés-Cabrera, under the GPLv2 and GPLv3 licences, together with the `libxdrfile` library of
GROMACS. That engine **was withdrawn from the project** and is not part of this distribution: no code, binary
or build dependency of cMMISMSA remains.

## External programs (not included)

QMeku invokes these programs if they are installed; it does not redistribute their code:

| Program | Use | Licence |
| --- | --- | --- |
| [xtb](https://github.com/grimme-lab/xtb) | GFN2-xTB, ALPB | LGPL-3.0 |
| [MOPAC](https://openmopac.net/) | PM6-D3H4X, COSMO | OpenMOPAC academic licence |
| [RDKit](https://www.rdkit.org/) | Ligands, hydrogens, formal charges | BSD-3-Clause |
| [PDB2PQR](https://github.com/Electrostatics/pdb2pqr) (optional) | Protonation of the receptor | BSD-3-Clause |

The scientific references of the methods (GFN2-xTB, PM6-D3H4X, ALPB, COSMO, MFCC and concaps) are in
`docs/METHOD.md` and in `README.md`.
