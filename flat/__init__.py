# Copyright (C) 2023-2025 Jure Cerar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
FLAT is a collection of scripts and utilities for PyMOL molecular visualizer.
"""

__author__ = "Jure Cerar"
__copyright__ = "Copyright (C) 2023-2025 Jure Cerar"
__license__ = "GNU GPL v3.0"
__version__ = "0.4.0"

try:
    from pymol import cmd
    pymol_version = cmd.get_version()[1]
except:
    pymol_version = "1.6"

__all__ = [
    "aaindex",
    "coloring",
    "creating",
    "editing",
    "electrostatics",
    "exporting",
    "fitting",
    "hydrophobic",
    "importing",
    "interface",
    "minimizing",
    "modelling",
    "plotting",
    "querying",
    "selecting",
    "seqalign",
    "shapes", 
    "viewing",
]

def init():
    """
    Imports all submodules and puts library into the pymol namespace
    for GUI menus. Also enables "help" in the PyMOL command line.
    """
    import sys
    import pymol
    from pymol import cmd

    # Init all submodules
    flat = __import__(__name__, fromlist=__all__)
    if not hasattr(pymol, "flat"):
        pymol.flat = flat

    # Allow import *
    if sys.modules.get("flat") != sys.modules[__name__]:
        sys.modules["flat"] = sys.modules[__name__]

    # pymol help
    if "flat" not in cmd.help_only:
        cmd.help_only["flat"] = [flat]
        cmd.help_sc.append("flat")

    return


def __init_plugin__(app=None):
    """ PyMOL Plugin hook """
    init()


# See: http://pymolwiki.org/index.php/Aa_codes
one_letter = {
    "PAQ": "Y", "AGM": "R", "ILE": "I", "PR3": "C", "GLN": "Q", "DVA": "V",
    "CCS": "C", "ACL": "R", "GLX": "Z", "GLY": "G", "GLZ": "G", "DTH": "T",
    "OAS": "S", "C6C": "C", "NEM": "H", "DLY": "K", "MIS": "S", "SMC": "C",
    "GLU": "E", "NEP": "H", "BCS": "C", "ASQ": "D", "ASP": "D", "SCY": "C",
    "SER": "S", "LYS": "K", "SAC": "S", "PRO": "P", "ASX": "B", "DGN": "Q",
    "DGL": "E", "MHS": "H", "ASB": "D", "ASA": "D", "NLE": "L", "DCY": "C",
    "ASK": "D", "GGL": "E", "STY": "Y", "SEL": "S", "CGU": "E", "ASN": "N",
    "ASL": "D", "LTR": "W", "DAR": "R", "VAL": "V", "CHG": "A", "TPO": "T",
    "CLE": "L", "GMA": "E", "HAC": "A", "AYA": "A", "THR": "T", "TIH": "A",
    "SVA": "S", "MVA": "V", "SAR": "G", "LYZ": "K", "BNN": "A", "5HP": "E",
    "IIL": "I", "SHR": "K", "HAR": "R", "FME": "M", "PYX": "C", "ALO": "T",
    "PHI": "F", "ALM": "A", "PHL": "F", "MEN": "N", "TPQ": "A", "GSC": "G",
    "PHE": "F", "ALA": "A", "MAA": "A", "MET": "M", "UNK": "X", "LEU": "L",
    "ALY": "K", "SET": "S", "GL3": "G", "TRG": "K", "CXM": "M", "TYR": "Y",
    "SCS": "C", "DIL": "I", "TYQ": "Y", "3AH": "H", "DPR": "P", "PRR": "A",
    "CME": "C", "IYR": "Y", "CY1": "C", "TYY": "Y", "HYP": "P", "DTY": "Y",
    "2AS": "D", "DTR": "W", "FLA": "A", "DPN": "F", "DIV": "V", "PCA": "E",
    "MSE": "M", "MSA": "G", "AIB": "A", "CYS": "C", "NLP": "L", "CYQ": "C",
    "HIS": "H", "DLE": "L", "CEA": "C", "DAL": "A", "LLP": "K", "DAH": "F",
    "HMR": "R", "TRO": "W", "HIC": "H", "CYG": "C", "BMT": "T", "DAS": "D",
    "TYB": "Y", "BUC": "C", "PEC": "C", "BUG": "L", "CYM": "C", "NLN": "L",
    "CY3": "C", "HIP": "H", "CSO": "C", "TPL": "W", "LYM": "K", "DHI": "H",
    "MLE": "L", "CSD": "A", "HPQ": "F", "MPQ": "G", "LLY": "K", "DHA": "A",
    "DSN": "S", "SOC": "C", "CSX": "C", "OMT": "M", "DSP": "D", "PTR": "Y",
    "TRP": "W", "CSW": "C", "EFC": "C", "CSP": "C", "CSS": "C", "SCH": "C",
    "OCS": "C", "NMC": "G", "SEP": "S", "BHD": "D", "KCX": "K", "SHC": "C",
    "C5C": "C", "HTR": "W", "ARG": "R", "TYS": "Y", "ARM": "R", "DNP": "A",
    "A": "A", "U": "U", "G": "G", "C": "C",  # RNA
    "DA": "A", "DT": "T", "DG": "G", "DC": "C",  # DNA
}

three_letter = {
    "A": "ALA", "C": "CYS", "E": "GLU", "D": "ASP", "G": "GLY", "F": "PHE",
    "I": "ILE", "H": "HIS", "K": "LYS", "M": "MET", "L": "LEU", "N": "ASN",
    "Q": "GLN", "P": "PRO", "S": "SER", "R": "ARG", "T": "THR", "W": "TRP",
    "V": "VAL", "Y": "TYR", "X": "UNK",
}
