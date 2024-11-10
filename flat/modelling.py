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
import collections
import sys
import io
from . import three_letter


@cmd.extend
def renumber(selection, start=1, quiet=False, *, _self=cmd):
    """
    DESCRIPTION
        Renumber sets new residue numbers (resi) for a polymer based on connectivity. 
    ARGUMENTS
        selection = string: atom selection to renumber {default: all}
        start = integer: counting start {default: 1}
    """
    # See: https://pymolwiki.org/index.php/Renumber
    start, quiet = int(start), bool(quiet)
    model = _self.get_model(selection)
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(10**5)
    _self.iterate(
        selection,
        "next(atom_it).model = model",
        space={
            "atom_it": iter(model.atom),
            "next": next,
        }
    )

    startatom = model.atom[0]
    for atom in model.atom:
        atom.adjacent = []
        atom.visited = False
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])
    minmax = [start, start]

    def traverse(atom, resi):
        atom.resi = resi
        atom.visited = True
        for other in atom.adjacent:
            if other.visited:
                continue
            if (atom.name, other.name) in [("C", "N"), ("O'", "P")]:
                minmax[1] = resi + 1
                traverse(other, resi + 1)
            elif (atom.name, other.name) in [("N", "C"), ("P", "O3'")]:
                minmax[0] = resi - 1
                traverse(other, resi - 1)
            elif (atom.name, other.name) not in [("SG", "SG")]:
                traverse(other, resi)
        return

    traverse(startatom, start)
    _self.alter(
        selection,
        "resi = next(atom_it).resi",
        space={
            "atom_it": iter(model.atom),
            "next": next,
        }
    )
    sys.setrecursionlimit(limit)
    if not quiet:
        print(" Renumber: range (%d to %d)" % tuple(minmax))
    return tuple(minmax)


@cmd.extend
def relax(selection, backbone=1, neighbors=1, cycles=100, state=0, *, _self=cmd):
    """
    DESCRIPTION
        Relax the given selection using sculpting wizard.
    USAGE
        relax selection [, backbone [, neighbors [, cycles [, state ]]]]
    """
    backbone, neighbors,  = int(backbone), int(neighbors)
    cycles, state = int(cycles), int(state)

    if not _self.get_object_list(selection):
        raise CmdException(f"Empty selection: {selection}")
    obj = _self.get_object_list(selection)[0]

    _self.protect()
    _self.deprotect(selection)
    if not backbone:
        _self.protect("name CA+C+N+O+OXT")
    if neighbors:
        _self.deprotect(f"byres ({obj} within 6.0 of ({selection}))")

    _self.sculpt_activate(obj, state)

    _self.set("sculpt_vdw_weight", 0.25, obj)  # Low VDW forces
    _self.set("sculpt_field_mask", 0x1FF, obj)  # Default
    _self.sculpt_iterate(obj, state, int(cycles * 0.75))

    _self.set("sculpt_field_mask", 0x01F, obj)  # Local Geometry Only
    _self.sculpt_iterate(obj, state, int(cycles * 0.25))

    _self.unset("sculpt_vdw_weight", obj)
    _self.unset("sculpt_field_mask", obj)
    _self.sculpt_deactivate(obj)
    _self.deprotect()

    return


@cmd.extend
def mutate(selection, residue, sculpt=0, cycles=100, *, _self=cmd):
    """
    DESCRIPTION
        Mutate a single residue and select best rotamer. Optionally
        fit lowest energy mutant rotamer to structure.
    USAGE
        mutate selection, residue [, sculpt [, cycles ]]
    """
    from . import three_letter
    sculpt, cycles = int(sculpt), int(cycles)
    
    if len(_self.get_model(selection).get_residues()) != 1:
        raise CmdException("Multiple residues in selection.")

    if len(residue) == 1:
        residue = three_letter.get(residue, residue)

    try:
        _self.wizard("mutagenesis")
        _self.do("refresh_wizard")
        _self.get_wizard().set_mode(residue.upper())
        _self.get_wizard().do_select(f"byres ({selection})")

        scores = _self.get_wizard().bump_scores
        if scores:
            frame = scores.index(min(scores)) + 1
        else:
            frame = 1
        _self.frame(str(frame))
        _self.get_wizard().apply()

    finally:
        _self.set_wizard()

    if sculpt > 0:
        relax(selection, 1, 1, cycles)

    return


@cmd.extend
def insert(selection, fragment, sculpt=0, cycles=100, *, _self=cmd):
    """
    DESCRIPTION
        Replace target residue with fragment. Use for non-standard residues. 
    USAGE
        insert target, frag [, out ]
    """
    sculpt, cycles = int(sculpt), int(cycles)

    if len(_self.get_model(selection).get_residues()) != 1:
        raise CmdException("Multiple residues in selection.")
    if len(_self.get_model(fragment).get_residues()) != 1:
        raise CmdException("Multiple residues in fragment.")

    temp = _self.get_unused_name("_tmp")
    _self.create(temp, fragment, source_state=1, target_state=1)

    model = _self.get_object_list(selection)[0]
    at = _self.get_model(selection).atom[0]
    segi, chain, resi = at.segi, at.chain, int(at.resi)
    _self.alter(temp, f"segi,chain,resi='{segi}','{chain}',{resi}")

    _self.align(
        mobile=f"({temp}) & n. O+C+N",
        target=f"({selection}) & n. O+C+N",
    )

    _self.remove(selection)
    _self.fuse(model, temp, mode=3)
    _self.edit()

    _self.bond(
        atom1=f"/{model}/{segi}/{chain}/`{resi-1}/C",
        atom2=f"/{model}/{segi}/{chain}/`{resi}/N",
    )
    _self.bond(
        atom1=f"/{model}/{segi}/{chain}/`{resi}/C",
        atom2=f"/{model}/{segi}/{chain}/`{resi+1}/N",
    )

    if sculpt > 0:
        relax(f"/{model}/{segi}/{chain}/`{resi}", 1, 1, cycles)

    return


@cmd.extend
def add_missing_atoms(selection, sculpt=0, cycles=100, *, _self=cmd):
    """
    DESCRIPTION
        Mutate those residues to themselves which have missing atoms.
    USAGE
        add_missing_atoms selection [, sculpt [, cycles ]]
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    """
    sculpt, cycles = int(sculpt), int(cycles)

    reference = {
        "ALA": {"CB"},
        "ARG": {"CB", "CG", "NE", "CZ", "NH1", "NH2", "CD"},
        "ASN": {"CB", "CG", "OD1", "ND2"},
        "ASP": {"CB", "CG", "OD1", "OD2"},
        "CYS": {"CB", "SG"},
        "GLN": {"CB", "CG", "CD", "NE2", "OE1"},
        "GLU": {"CB", "CG", "OE2", "CD", "OE1"},
        "GLY": set(),
        "HIS": {"CE1", "CB", "CG", "CD2", "ND1", "NE2"},
        "ILE": {"CB", "CD1", "CG1", "CG2"},
        "LEU": {"CB", "CG", "CD1", "CD2"},
        "LYS": {"CB", "CG", "NZ", "CE", "CD"},
        "MET": {"CB", "CG", "CE", "SD"},
        "PHE": {"CE1", "CB", "CG", "CZ", "CD1", "CD2", "CE2"},
        "PRO": {"CB", "CG", "CD"},
        "SER": {"OG", "CB"},
        "THR": {"CB", "OG1", "CG2"},
        "TRP": {"CZ2", "CB", "CG", "CH2", "CE3", "CD1", "CD2", "CZ3", "NE1", "CE2"},
        "TYR": {"CE1", "OH", "CB", "CG", "CZ", "CD1", "CD2", "CE2"},
        "VAL": {"CB", "CG1", "CG2"},
    }

    namelists = collections.defaultdict(list)
    _self.iterate(
        f"({selection}) & polymer",
        "namelists[model,segi,chain,resn,resi].append(name)",
        space=locals()
    )

    sele_list = []
    for key, namelist in namelists.items():
        resn = key[3]
        if resn not in reference:
            print(f"Unknown residue: '{resn}'")
            continue

        if not reference[resn].issubset(namelist):
            try:
                _self.wizard("mutagenesis")
                _self.do("refresh_wizard")
                _self.get_wizard().set_mode(resn)
                sele = "/{}/{}/{}/{}`{}".format(*key)
                _self.get_wizard().do_select(sele)
                scores = _self.get_wizard().bump_scores
                if scores:
                    frame = scores.index(min(scores)) + 1
                else:
                    frame = 1
                _self.frame(str(frame))
                _self.get_wizard().apply()
                sele_list.append(sele)
                # if not quiet:
                #     print(f"Mutated: {sele}")

            except Exception as ex:
                print(f"Mutating '{sele}' failed: {ex}")

            finally:
                _self.set_wizard()

    if sculpt:
        for sele in sele_list:
            relax(sele, 0, 0, cycles)

    return 


@cmd.extend
def fix_h(selection="all", *, _self=cmd):
    """
    DESCRIPTION
        Adds hydrogens with correct names onto a molecule using OpenMM based on known templates.
    USAGE
        fix_h [selection]
    PARAMETERS
        selection = str: atom selection {default: all}
    NOTES
        It's a kludge and may need some love and persuasion to work. This adds hydrogens based
        on template i.e. it will most likely only work for proteins, sugars, and DNA/RNA. Also
        needs to have well defined C- and N-terminus.  
    """
    import openmm
    import openmm.app

    if _self.get_model(selection).nAtom == 0:
        raise CmdException("Empty selection")

    objects = _self.get_object_list(selection)
    if len(objects) != 1:
        raise CmdException("Multiple objects in selection")

    # Pass structure to OpenMM with no hydrogens
    molstr = _self.get_str("pdb", f"({selection}) and ! hydrogen", -1)
    with io.StringIO(molstr) as f:
        structure = openmm.app.PDBFile(f)

    # Use modeller to add hydrogens from template
    modeller = openmm.app.Modeller(structure.topology, structure.positions)
    modeller.addHydrogens()

    # Load structure back
    with io.StringIO() as f:
        openmm.app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        molstr = f.getvalue()

    # Replace original structure
    tmp = _self.get_unused_name("_tmp")
    try:
        _self.load_raw(molstr, "pdb", tmp, 0, zoom=0)
        _self.remove(selection)
        _self.copy_to(objects[0], tmp, "", zoom=0)
    finally:
        _self.delete(tmp)


@cmd.extend
def update_align(mobile, target, state=1, *, fix="none", quiet=1, _self=cmd):
    """
    DESCRIPTION
        Update (and optionally fix) coordinates based on sequence alignment.
    USAGE
        update_align mobile, target [, state [, fix ]]
    PARAMETERS
        mobile = str: atom selection of mobile section.
        target = str: atom selection of target section. 
        state = int: select state. {default: 1}
        fix = restrain | fix | protect | none: Method for fixing updated atoms
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    """
    aln = _self.get_unused_name("aln_hom")
    _self.align(mobile, target, object=aln, cycles=0, max_gap=-1)
    try:
        mobile_aln = f"({mobile}) & {aln}"
        _self.update(mobile_aln,
                     f"({target}) & {aln}",
                     state,
                     state,
                     matchmaker=0,
                     quiet=quiet)
        if fix == "restrain":
            _self.reference("store", mobile_aln, state, quiet=quiet)
            _self.flag("restrain", mobile_aln, "set", quiet=quiet)
        elif fix == "fix":
            _self.flag("fix", mobile_aln, "set", quiet=quiet)
        elif fix == "protect":
            _self.protect(mobile_aln, quiet=quiet)
        elif fix != "none":
            raise ValueError(fix)
    finally:
        _self.delete(aln)


@cmd.extend
def sculpt_homolog(mobile, target, state=1, cycles=1000, *, fix="restrain", quiet=1, _self=cmd):
    """
    DESCRIPTION
        Sculpt mobile towards target, based on sequence alignment.
    USAGE 
        sculpt_homolog mobile, target [, state [, cycles [, fix ]]]
    ARGUMENTS
        mobile = str: atom selection of mobile section.
        target = str: atom selection of target section. 
        state = int: select state. {default: 1}
        cycles: Number of sculpt iterations. {default: 1000}
        fix = restrain | fix | protect | none: Method for fixing updated atoms
    SOURCE
        From PSICO(c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    """
    (mobile_object, ) = _self.get_object_list(mobile)
    _self.sculpt_activate(mobile_object, state)
    update_align(mobile, target, state, fix=fix, quiet=quiet, _self=_self)
    _self.sculpt_iterate(mobile, state, cycles)


# Autocomplete
cmd.auto_arg[0].update({
    "mutate": cmd.auto_arg[0]["zoom"],
    "relax": cmd.auto_arg[0]["zoom"],
    "mutate": cmd.auto_arg[0]["zoom"],
    "insert": cmd.auto_arg[0]["zoom"],
    "add_missing_atoms": cmd.auto_arg[0]["zoom"],
    "fix_h": cmd.auto_arg[0]["zoom"],
    "update_align": cmd.auto_arg[0]["align"],
    "sculpt_homolog": cmd.auto_arg[0]["align"],
})
cmd.auto_arg[1].update({
    "mutate": [cmd.Shortcut(three_letter.values()), "residue", ""],
    "update_align": cmd.auto_arg[1]["align"],
    "sculpt_homolog": cmd.auto_arg[1]["align"],
})

