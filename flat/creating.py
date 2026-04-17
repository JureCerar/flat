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
:mod:`flat.creating`
====================
Module for creating new molecular objects.
"""

from pymol import cmd, CmdException
from . import three_letter

# Cheat sheet for atom side chain fixing 
_SIDECHAIN_CENTER_ATOMS = {
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

# Available side chain center methods
_SIDECHAIN_CENTER_METHODS = ["bahar1996", "centroid"]


# Atoms for determining center of aromatic residues
_AROMATIC_CENTER_ATOMS = { 
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

        1.  Sidechain interaction centers as defined by Bahar and Jernigan 1996
            http://www.ncbi.nlm.nih.gov/pubmed/9080182  
        2.  Sidechain centroids, the pseudoatom is the centroid of all atoms except
            hydrogens and backbone atoms (N, C and O).  

        With method `bahar1996`, if a residue has all relevant side-chain center
        atoms missing (for example a MET without SD), it will be missing in the
        created pseudo-atom object.

        With method `centroid`, if you want to exclude C-alpha atoms from
        side-shains.

    USAGE
        sidechain_centers [ selection [, object [, method, [ name ]]]]
    ARGUMENTS
        selection : str, optional
            Atom selection. 
        object : str, default = 'sidechain_centers'
            Name of object to create.
        method : str, default = 'bahar1996'
            Method for calculating residue centroid.
        name : str, default = 'PS1'
            Atom name of created pseudo-atoms.
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder
    """
    from pymol.chempy import Atom, cpv, models

    atmap = dict()
    if method in ["bahar1996", "1", 1]:
        modelAll = _self.get_model("(%s) and resn %s" % (
            selection, "+".join(_SIDECHAIN_CENTER_ATOMS)))
        for at in modelAll.atom:
            if at.name in _SIDECHAIN_CENTER_ATOMS[at.resn]:
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
        Creates an object with pseudo-atoms representing aromatic centers for 
        each residue in selection.
    USAGE
        aromatic_centers [selection [, object [, name ]]]
    ARGUMENTS
        selection : str, optional
            Atom selection. 
        object : str, default = 'aromatic_centers'
            Name of object to create.
        name : str, default = 'PS1'
            Atom name of created pseudo-atoms.
    RETURNS
        : chempy.model
            Returns molecular model for aromatic centers.
    """
    from pymol.chempy import Atom, cpv, models

    atmap = dict()

    res_names = "+".join(_AROMATIC_CENTER_ATOMS)
    modelAll = _self.get_model(f"({selection}) and resn {res_names}")
    for at in modelAll.atom:
        if at.name in _AROMATIC_CENTER_ATOMS[at.resn]:
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
def com(selection="all", name="COM", *, _self=cmd):
    """
    DESCRIPTION
        Calculates the center of mass. Considers atom mass and occupancy.
    USAGE
        com [ selection [, name ]]
    ARGUMENTS
        selection : str, optional
            Atom selection. 
        name : str, default = 'COM'
            Name of object to create.
    """
    if not name:
        name = _self.get_unused_name("COM")
    for i in range(_self.count_states(selection)):
        x, y, z = _self.centerofmass(selection, state=i+1)
        _self.pseudoatom(name, pos=[x, y, z], state=i+1)
    _self.show("nonbonded", name)


@cmd.extend
def fragment(name, object=None, origin=1, zoom=0, quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Retrieves a 3D structure from the fragment library, which is 
        currently pretty meager.

        This overload internal `fragment` function.
    USAGE
        fragment name [, object [, origin [, zoom ]]]
    ARGUMENTS
        name : str
            Name of library fragment
        object : str, default = None
            Name of object to create. Default is fragment name.
        origin : bool, default = True
            Center fragment at the current position.
        zoom : bool, default = False
            Zoom view to fit fragment
    RETURNS
        : chempy.model
            Chempy model of loaded fragment. 
    """
    # NOTE: Do not change function signature
    import chempy
    import os

    if object is None:
        object = name

    r = cmd.DEFAULT_ERROR

    # First try to get model from default library
    chempy_path = os.path.join(chempy.path, "fragments", name + ".pkl")
    flatcat_path = os.path.join(os.path.dirname(__file__),
                                "fragments", name + ".pkl")
    if os.path.isfile(chempy_path):
        path = chempy_path
    elif os.path.isfile(flatcat_path):
        path = flatcat_path
    else:
        raise CmdException(f"Unable to load fragment: '{name}'")

    model = chempy.io.pkl.fromFile(path)
    la = len(model.atom)
    if la and int(origin):
        position = _self.get_position()
        for c in range(0, 3):
            mean_c = sum([a.coord[c] for a in model.atom]) / la
            mean_c = position[c] - mean_c
            for a in model.atom:
                a.coord[c] += mean_c
    r = _self.load_model(model, str(object), quiet=quiet,
                         zoom=zoom, _self=_self)
    return r


# Autocomplete
cmd.auto_arg[0].update({
    "sidechain_centers": cmd.auto_arg[0]["zoom"],
    "aromatic_centers": cmd.auto_arg[0]["zoom"],
    "com": cmd.auto_arg[0]["zoom"],
    "fragment": [cmd.Shortcut(three_letter.values()), "residue", ""],
})

cmd.auto_arg[2].update({
    "sidechain_centers": [cmd.Shortcut(_SIDECHAIN_CENTER_METHODS), "method", ""],
})
