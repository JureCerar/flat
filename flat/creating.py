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

sidechain_center_atoms = {
    "GLY": ("CA",),
    "ALA": ("CB",),
    "VAL": ("CG1", "CG2"),
    "ILE": ("CD1",),
    "LEU": ("CD1", "CD2"),
    "SER": ("OG",),
    "THR": ("OG1", "CG2"),
    "ASP": ("OD1", "OD2"),
    "ASN": ("OD1", "ND2"),
    "GLU": ("OE1", "OE2"),
    "GLN": ("OE1", "NE2"),
    "LYS": ("NZ",),
    "ARG": ("NE", "NH1", "NH2"),
    "CYS": ("SG",),
    "MET": ("SD",),
    "MSE": ("SE",),
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "PRO": ("CB", "CG", "CD"),
}

sidechain_center_methods = ["bahar1996", "centroid"]

aromatic_center_atoms = { 
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ",),
}


@cmd.extend
def sidechain_centers(selection="all", object="sidechain_centers", method="bahar1996", name="PS1", *, _self=cmd):
    """
    DESCRIPTION
        Creates an object with sidechain representing pseudoatoms for each residue
        in selection.

        Two methods are available:
        (1) Sidechain interaction centers as defined by Bahar and Jernigan 1996
            http://www.ncbi.nlm.nih.gov/pubmed/9080182
        (2) Sidechain centroids, the pseudoatom is the centroid of all atoms except
            hydrogens and backbone atoms (N, C and O).

        With method "bahar1996", if a residue has all relevant sidechain center
        atoms missing (for example a MET without SD), it will be missing in the
        created pseudoatom object.

        With method "centroid", if you want to exclude C-alpha atoms from
        sidechains, modify the selection like in this example:

        >>> sidechain_centers newobject, all and (not name CA or resn GLY), method=2
    USAGE
        sidechain_centers [object [, selection [, method, [ name ]]]]
    ARGUMENTS
        object = string: name of object to create
        selection = string: atoms to consider {default: (all)}
        method = string: bahar1996 or centroid {default: bahar1996}
        name = string: atom name of pseudoatoms {default: PS1}
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder
    """
    from pymol.chempy import Atom, cpv, models

    atmap = dict()
    if method in ["bahar1996", "1", 1]:
        modelAll = _self.get_model("(%s) and resn %s" % (
            selection, "+".join(sidechain_center_atoms)))
        for at in modelAll.atom:
            if at.name in sidechain_center_atoms[at.resn]:
                atmap.setdefault(
                    (at.segi, at.chain, at.resn, at.resi), []).append(at)
    elif method in ["centroid", "2", 2]:
        modelAll = _self.get_model(
            "(%s) and polymer and not (hydro or name C+N+O)" % selection)
        for at in modelAll.atom:
            atmap.setdefault(
                (at.segi, at.chain, at.resn, at.resi), []).append(at)
    else:
        raise CmdException("unknown method: {}".format(method))

    model = models.Indexed()
    for centeratoms in atmap.values():
        center = cpv.get_null()
        for at in centeratoms:
            center = cpv.add(center, at.coord)
        center = cpv.scale(center, 1. / len(centeratoms))
        atom = Atom()
        atom.coord = center
        atom.index = model.nAtom + 1
        atom.name = name
        for key in ["segi", "chain", "resi_number", "resi", "resn", "hetatm", "ss", "b",]:
            setattr(atom, key, getattr(at, key))
        model.add_atom(atom)
    model.update_index()
    if object in _self.get_object_list():
        _self.delete(object)
    _self.load_model(model, object)
    return model


@cmd.extend
def aromatic_centers(selection="all", object="aromatic_centers", name="PS1", *, _self=cmd):
    """
    DESCRIPTION
        Creates an object with pseudoatoms representing atomatic centers for 
        each residue in selection.
    USAGE
        aromatic_centers [object [, selection [, name ]]]
    ARGUMENTS
        object = string: name of object to create
        selection = string: atoms to consider {default: (all)}
        name = string: atom name of pseudoatoms {default: PS1}
    """
    from pymol.chempy import Atom, cpv, models

    atmap = dict()

    res_names = "+".join(aromatic_center_atoms)
    modelAll = _self.get_model(f"({selection}) and resn {res_names}")
    for at in modelAll.atom:
        if at.name in aromatic_center_atoms[at.resn]:
            atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)

    model = models.Indexed()
    for aromatic_atoms in atmap.values():
        center = cpv.get_null()
        for at in aromatic_atoms:
            center = cpv.add(center, at.coord)
        center = cpv.scale(center, 1. / len(aromatic_atoms))
        atom = Atom()
        atom.coord = center
        atom.index = model.nAtom + 1
        atom.name = name
        for key in ["segi", "chain", "resi_number", "resi", "resn", "hetatm", "ss", "b",]:
            setattr(atom, key, getattr(at, key))
        model.add_atom(atom)
    model.update_index()
    if object in _self.get_object_list():
        _self.delete(object)
    _self.load_model(model, object)
    return model


@cmd.extend
def com(selection="all", object=None, state=0, *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Calculates the center of mass. Considers atom mass and occupancy.
    USAGE
        com [object [, selection [, name ]]]
    ARGUMENTS
        selection = string: atoms to consider {default: (all)}
        object = string: name of object to create {default: None}
        state = int: object state, -1 for current state, 0 for all states {default: 0}
    """
    state, quiet = int(state), int(quiet)
    if (object == None):
        try:
            object = _self.get_legal_name(selection)
            object = _self.get_unused_name(object + "_COM", 0)
        except AttributeError:
            object = "COM"
    _self.delete(object)

    if (state != 0):
        x, y, z = _self.centerofmass(selection, state=state)
        _self.pseudoatom(object, pos=[x, y, z])
        _self.show("spheres", object)
        
    else:
        for i in range(_self.count_states()):
            x, y, z = _self.centerofmass(selection, state=i+1)
            _self.pseudoatom(object, pos=[x, y, z], state=i+1)
            _self.show("spheres", object)


# Autocomplete
cmd.auto_arg[0].update({
    "sidechain_centers": cmd.auto_arg[0]["zoom"],
    "aromatic_centers": cmd.auto_arg[0]["zoom"],
    "com": cmd.auto_arg[0]["zoom"],
})

cmd.auto_arg[2].update({
    "sidechain_centers": [cmd.Shortcut(sidechain_center_methods), "method", ""],
})
