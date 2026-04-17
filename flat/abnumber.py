# Copyright (C) 2023-2026 Jure Cerar
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
:mod:`flat.abnumber`
====================
Module for Antibody numbering and alignment using ANARCI

Installation
------------
This is an optional dependency of flat and can be installed with:

.. code:: bash

    pip install abnumber

.. note::

    Windows is not supported due to `HMMER`_ dependency. AbNumber
    is currently only available on UNIX & MacOS.

.. _HMMER:
    https://github.com/EddyRivasLab/hmmer
"""


from pymol import cmd, CmdException

# Antibody numbering schemes
_CDR_SCHEMES = ["imgt", "chothia", "kabat", "aho"]
"""Antibody numbering schemes"""
 

@cmd.extend
def abnumber(selection="all", scheme="chothia", *, _self=cmd):
    """
    DESCRIPTION
        Select antibody CDR region according to Antibody Numbering
        and Antigen Receptor Classification (ANARCI).
    USAGE
        abnumber [ selection [, scheme ]]
    ARGUMENTS
        selection : str, optional
            Atom selection.
        scheme : str, default = 'chothia'
            Which numbering scheme should be used.
    SCHEMES
        imgt, kabat, chothia, aho
    SEE ALSO
        https://github.com/prihoda/AbNumber
    """
    from abnumber import Chain
    from abnumber.exceptions import ChainParseError

    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple objects in selection")
    if scheme not in _CDR_SCHEMES:
        raise CmdException(f"Invalid numbering scheme: '{scheme}'")

    for chain in _self.get_chains(selection):
        subsele = f"({selection}) & c. {chain}"

        # Get sequence
        fasta = cmd.get_fastastr(subsele)
        seq = "".join(fasta.splitlines()[1:])

        try:
            cdr = Chain(seq, scheme, use_anarcii=True)
        except ChainParseError:
            print("Warning: Variable chain sequence not recognized")
            continue

        # Create selection
        name = _self.get_unused_name(f"{chain}_CDR", 0)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr1_seq}", 0, merge=1)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr2_seq}", 0, merge=1)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr3_seq}", 0, merge=1)


@cmd.extend
def abrenumber(selection="(all)", scheme="chothia", *, _self=cmd):
    """
    DESCRIPTION
        Renumber and rename (chains) antibodies according to 
        Antibody Numbering and Antigen Receptor Classification (ANARCI)
    USAGE
        abrenumber [ selection [, scheme ]]
    ARGUMENTS
        selection : str, optional
            Atom selection.
        scheme : str, default = 'chothia'
            Which numbering scheme should be used.
    SCHEMES
        imgt, kabat, chothia, aho
    SEE ALSO
        https://github.com/prihoda/AbNumber
    """
    import sys
    from abnumber import Chain
    from abnumber.exceptions import ChainParseError

    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple objects in selection")
    if scheme not in _CDR_SCHEMES:
        raise CmdException(f"Invalid numbering scheme: '{scheme}'")

    # Raise recursion limit for large molecules
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(10**5)

    hc_name, lc_name = "H", "L"

    for chain in _self.get_chains(selection):
        subsele = f"({selection}) & c. {chain}"

        # Get sequence
        fasta = _self.get_fastastr(subsele)
        seq = "".join(fasta.splitlines()[1:])

        try:
            cdr = Chain(seq, scheme, use_anarcii=True)
        except ChainParseError:
            print("Warning: Variable chain sequence not recognized")
            continue

        # Get ordered list of residue indices
        abnum = []
        for pos, aa in cdr:
            abnum.append(str(pos)[1:])
        last = int(abnum[-1])

        # Alter model for BFS
        # See: https://pymolwiki.org/index.php/Renumber
        model = _self.get_model(subsele)
        startatom = model.atom[0]
        for atom in model.atom:
            atom.adjacent = []
            atom.visited = False
        for bond in model.bond:
            atoms = [model.atom[i] for i in bond.index]
            atoms[0].adjacent.append(atoms[1])
            atoms[1].adjacent.append(atoms[0])

        # Traverse selection and increment index if we pass C-N bond
        # between two peptides.
        def traverse(atom, i):
            if i < len(abnum):
                atom.resi = abnum[i]
            else:
                atom.resi = last + (i - len(abnum))
            atom.visited = True
            for other in atom.adjacent:
                if other.visited:
                    continue
                if (atom.name, other.name) in [("C", "N"), ("O'", "P")]:
                    traverse(other, i + 1)
                elif (atom.name, other.name) in [("N", "C"), ("P", "O3'")]:
                    traverse(other, i - 1)
                elif (atom.name, other.name) not in [("SG", "SG")]:
                    traverse(other, i)
            return

        # Let's do this
        traverse(startatom, 0)
        _self.alter(
            subsele,
            "resi = next(atom_it).resi",
            space={
                "atom_it": iter(model.atom),
                "next": next,
            }
        )

        # Alter chain names
        if cdr.is_heavy_chain():
            _self.alter(subsele, f"chain='{hc_name}'")
            hc_name = chr(ord(hc_name) + 1)
        else:
            _self.alter(subsele, f"chain='{lc_name}'")
            lc_name = chr(ord(lc_name) + 1)

    # Set back original recursion limit
    # that we don't break anything
    sys.setrecursionlimit(limit)


# Autocomplete
cmd.auto_arg[0].update({
    "abnumber": cmd.auto_arg[0]["zoom"],
    "abrenumber": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "abnumber": [cmd.Shortcut(_CDR_SCHEMES), "scheme", ""],
    "abrenumber": [cmd.Shortcut(_CDR_SCHEMES), "scheme", ""],
})
