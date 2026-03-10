# Copyright (C) 2023-2024 Jure Cerar
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


@cmd.extend
def find_seq(pattern, selection="all", name="sele", merge=False, *, _self=cmd):
    """
    DESCRIPTION
        Given a sequence/regex to find, select those matching 
        amino acids in the protein.
    USAGE
        find_seq pattern, [ selection [, name [, merge ]]]
    ARGUMENTS
        pattern = str: the sequence of amino acids to match and selection.
            This can be a sequence of amino acids or a regular expression.  
        selection = str: a selection-expression. {default: "all"}
        name = str: a unique name for the selection. {default: "sele"}
        merge = int: merge to existing selection. {default: False}
    EXAMPLE
        >>> find_seq N[^P][TS]
    """
    from . import one_letter
    import re

    merge = bool(merge)

    if len(pattern) == 0 or not isinstance(pattern, str):
        raise CmdException("Search sequence not provided")
    if len(selection) == 0 or not isinstance(selection, str):
        raise CmdException("Invalid selection or object")
    if len(name) == 0 or not isinstance(name, str):
        raise CmdException("Invalid name for the selection")
    
    regex = re.compile(pattern.upper())

    # Merge or create new selection
    _self.select(name, "None", merge=merge)

    for object in _self.get_object_list(selection):
        for chain in _self.get_chains(f"o. {object}"):
            residues = []
            _self.iterate(
                f"bca. ({selection}) and o. {object} and c. {chain}",
                "residues.append((resi, resn))",
                space=locals()
            )
            sequence = "".join([one_letter.get(r[1], "-") for r in residues])
            for m in regex.finditer(sequence):
                (start, stop) = m.span()
                indices = "+".join(residues[i][0] for i in range(start, stop))
                _self.select(
                    name,
                    f"({selection}) and o. {object} and c. {chain} and i. {indices}",
                    merge=True
                )

    return name


@cmd.extend
def diff(sele1, sele2, byres=1, name=None, operator="in", quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Difference between two molecules
    USAGE 
        diff sele1, sele2 [, byres [, name [, operator ]]]
    ARGUMENTS
        sele1 = string: atom selection.
        sele2 = string: atom selection.
        byres = 0/1: report residues, not atoms (does not affect selection). {default: 1}
        operator = in/like/align: operator to match atoms {default: in}
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    SEE ALSO
        symdiff
    """
    byres, quiet = int(byres), int(quiet)
    if name is None:
        name = _self.get_unused_name("diff")
    if operator == "align":
        alnobj = _self.get_unused_name("__aln")
        _self.align(sele1, sele2, cycles=0, transform=0, object=alnobj)
        sele = "(%s) and not %s" % (sele1, alnobj)
        _self.select(name, sele)
        _self.delete(alnobj)
    else:
        sele = "(%s) and not ((%s) %s (%s))" % (sele1, sele1, operator, sele2)
        _self.select(name, sele)
    if not quiet:
        if byres:
            seleiter = "byca " + name
            expr = "print('/%s/%s/%s/%s`%s' % (model,segi,chain,resn,resi))"
        else:
            seleiter = name
            expr = "print('/%s/%s/%s/%s`%s/%s' % (model,segi,chain,resn,resi,name))"
        _self.iterate(seleiter, expr)
    return name


@cmd.extend
def symdiff(sele1, sele2, byres=1, name=None, operator="in", quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Symmetric difference between two molecules
    USAGE 
        symdiff sele1, sele2 [, byres [, name [, operator ]]]
    ARGUMENTS
        sele1 = string: atom selection.
        sele2 = string: atom selection.
        byres = 0/1: report residues, not atoms (does not affect selection). {default: 1}
        operator = in/like/align: operator to match atoms {default: in}
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    SEE ALSO
        diff
    """
    byres, quiet = int(byres), int(quiet)
    if name is None:
        name = _self.get_unused_name("symdiff")
    tmpname = _self.get_unused_name("__tmp")
    diff(sele1, sele2, byres, name, operator, quiet)
    diff(sele2, sele1, byres, tmpname, operator, quiet)
    _self.select(name, tmpname, merge=1)
    _self.delete(tmpname)
    return name


def wait_for(name, state=0, quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Wait for "name" to be available as selectable object.
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    """
    if _self.count_atoms("?" + name, 1, state) == 0:
        s = _self.get_setting_boolean("suspend_updates")
        if s:
            _self.set("suspend_updates", 0)
        _self.refresh()
        if s:
            _self.set("suspend_updates")


# Autocomplete
cmd.auto_arg[0].update({
    "diff": cmd.auto_arg[0]["align"],
    "symdiff": cmd.auto_arg[0]["align"],
})
cmd.auto_arg[1].update({
    "find_seq": cmd.auto_arg[1]["select"],
    "diff": cmd.auto_arg[1]["align"],
    "symdiff": cmd.auto_arg[1]["align"],
})
