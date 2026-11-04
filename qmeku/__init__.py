# QMeku: approximate binding energy and per-residue decomposition of docking poses.
# Copyright (C) 2026 SB2CS lab
# SPDX-License-Identifier: GPL-3.0-or-later
# Free software under the GNU General Public License, either version 3 or (at
# your option) any later version; see the LICENSE file.
# Distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
"""QMeku: approximate binding energy and per-residue decomposition of docking poses with semiempirical methods
(xtb, MOPAC), including the gas / solvation balance and the GFN2 energy components."""
__version__ = "1.0"

from .analysis import InputError, analyze_ligand, analyze_pose, pocket_residues  # noqa: F401
from .ligand import QmekuError, build_complex, prepare_ligand, prepare_receptor, read_poses  # noqa: F401
from .report import format_log, summarize_residues, write_outputs  # noqa: F401
