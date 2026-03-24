from pymol import cmd
from pymol import editor

import flat.creating

aa_dict = {
    "A": "ALA", "C": "CYS", "E": "GLU", "D": "ASP", "G": "GLY", "F": "PHE",
    "I": "ILE", "H": "HIS", "K": "LYS", "M": "MET", "L": "LEU", "N": "ASN",
    "Q": "GLN", "P": "PRO", "S": "SER", "R": "ARG", "T": "THR", "W": "TRP",
    "V": "VAL", "Y": "TYR",
    # Custom residues
    "x": "OCS", "y": "CSD"
}

# IMPORTANT: Overload pymol build in with alternative
cmd.fragment = flat.creating.fragment

# Object name and sequence
name = "obj01"
sequence = "GGxGGyGG"
ss = 4  # 1=alpha helix, 2=antiparallel beta, 3=parallel beta, 4=flat

# Build peptide
code = sequence[0]
cmd.fragment(aa_dict[code], name)
cmd.alter(name, "resi = 1")
cmd.edit(f"{name} and name C")
for code in sequence[1:]:
    editor.attach_amino_acid("pk1", aa_dict[code], ss=ss)
cmd.edit()
cmd.zoom()
