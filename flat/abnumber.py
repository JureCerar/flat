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

from pymol import cmd, CmdException


_CDR_SCHEMES = ["imgt", "chothia", "kabat", "aho"]
"""Antibody numbering schemes"""


@cmd.extend
def abnumber(selection="(all)", scheme="imgt", *, _self=cmd):
    """
    DESCRIPTION
        Antibody numbering and antigen receptor classification
    USAGE
        abnumber [ selection [, scheme [, exe ]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        scheme = str: Which numbering scheme should be used. {default: imgt}
    SCHEMES
        imgt, kabat, chothia, aho
    SEE ALSO
        ...
    """
    from . import one_letter
    from abnumber import Chain
    from abnumber.exceptions import ChainParseError

    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple objects in selection")
    if scheme not in _CDR_SCHEMES:
        raise CmdException(f"Invalid numbering scheme: '{scheme}'")

    def seqbuilder(seq_list):
        result = ""
        for resn in seq_list:
            if resn in one_letter:
                result += one_letter[resn]
            else:
                print(f"Warning: Unknown residue '{resn}'")
                result += "-"
        return result

    for chain in _self.get_chains(selection):
        seq_list = []
        _self.iterate(
            f"(bca. {selection} & c. {chain}) & polymer",
            "seq_list.append(resn)",
            space=locals(),
        )

        try:
            cdr = Chain(seqbuilder(seq_list), scheme)
        except ChainParseError:
            print("Warning: Variable chain sequence not recognized")
            continue

        subsele = f"{selection} & c. {chain}"
        name = _self.get_unused_name(f"{chain}_CDR", 0)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr1_seq}", 0, merge=1)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr2_seq}", 0, merge=1)
        _self.select(name, f"{subsele} & pepseq {cdr.cdr3_seq}", 0, merge=1)


# Autocomplete
cmd.auto_arg[0].update({
    "abnumber": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "abnumber": [cmd.Shortcut(_CDR_SCHEMES), "scheme", ""],
})
