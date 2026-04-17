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
:mod:`flat.editing`
===================
Module for editing molecular objects.
"""

from pymol import cmd, CmdException
import numpy as np
import sys

# Available options for atom unwrapping
_UNWRAP_MODES = ["residues", "chains", "segments", "bonds"]


@cmd.extend
def renumber(selection, start=1, *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Renumber sets new residue numbers (`resi`) for 
        a polymer based on connectivity.
    USAGE
        renumber selection [, start ]
    ARGUMENTS
        selection : str
            Atom selection.
        start : int, default = 1
            Numbering start.
    RETURNS
        : Tuple(int, int)
            First and last index in renumbering.  
    """
    # See: https://pymolwiki.org/index.php/Renumber
    start, quiet = int(start), bool(quiet)
    model = _self.get_model(selection)
    if not model.atom:
        raise CmdException("Invalid selection")

    # NOTE: Increase recursion limit for large molecules
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

    # Traverse model and increase residue number if
    # we pass specific bonds (peptide, disulfide, etc.)
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

    # Let's do this
    traverse(startatom, start)
    _self.alter(
        selection,
        "resi = next(atom_it).resi",
        space={
            "atom_it": iter(model.atom),
            "next": next,
        }
    )

    # Set back original recursion limit
    sys.setrecursionlimit(limit)

    if not quiet:
        print("Renumber: range (%d to %d)" % tuple(minmax))

    return tuple(minmax)


@cmd.extend
def wrap(selection="all", *, _self=cmd):
    """
    DESCRIPTION
        Shift the contents of a given selection back into the unit cell::

            +-----------+            +-----------+
            |       1-2-|-3-4  -->   |-3-4   1-2-|
            +-----------+            +-----------+

    USAGE
        wrap [ selection ]
    ARGUMENTS
        selection : str, optional
            Atoms selection. 
    """
    states = _self.count_states(selection)
    for state in range(1, states + 1):
        box = _self.get_symmetry(selection, state)[0:3]
        xyz = _self.get_coords(selection, state)
        # Apply PBC and update coordinates
        _self.load_coords(xyz % box, selection, state)


@cmd.extend
def unwrap(selection="all", mode="chains", *, _self=cmd):
    """
    DESCRIPTION
        Move all atoms in an selection so that bonds don't split over images::

            +-----------+        +-----------+
            |-3-4   1-2-|  -->   |       1-2-|-3-4
            +-----------+        +-----------+

    USAGE
        unwrap [ selection [, mode ]]
    ARGUMENTS
        selection : str, optional
            Atoms selection. 
        mode : str, default = 'chains'
            The group which will be kept together through the shifting
            process: `residues`, `chains`, `segments`, or `bonds`.
    """
    # Get model object and modify it
    model = _self.get_model(selection)
    for i, atom in enumerate(model.atom):
        # Ensure atom index matches it's array index
        atom.index = i
        atom.adjacent = []
    for bond in model.bond:
        atoms = [model.atom[i] for i in bond.index]
        atoms[0].adjacent.append(atoms[1])
        atoms[1].adjacent.append(atoms[0])

    if mode == "bonds":
        # Get list of atoms indices that are bonded together
        VISITED = set()  # Atoms we visited

        def traverse(atom):
            """Traverse bonded network"""
            VISITED.add(atom.index)
            network = {atom.index}
            for other in atom.adjacent:
                if other.index in VISITED:
                    continue
                nt = traverse(other)
                network.update(nt)
            return network

        try:
            # Increase recursion limit
            limit = sys.getrecursionlimit()
            sys.setrecursionlimit(10**5)
            network = []
            for atom in model.atom:
                if atom.index in VISITED:
                    continue
                nt = traverse(atom)
                network.append(list(nt))
        finally:
            sys.setrecursionlimit(limit)

    elif mode in ["res", "residues"]:
        # Get list of atoms indices that belong to same residue
        residues = dict()
        for atom in model.atom:
            key = (atom.segi, atom.chain, atom.resi)
            if key in residues:
                residues[key].append(atom.index)
            else:
                residues[key] = [atom.index]
        network = residues.values()

    elif mode in ["c", "chains"]:
        # Get list of atoms indices that belong to same chain
        chains = dict()
        for atom in model.atom:
            key = (atom.segi, atom.chain)
            if key in chains:
                chains[key].append(atom.index)
            else:
                chains[key] = [atom.index]
        network = chains.values()

    elif mode in ["seg", "segments"]:
        # Get list of atoms indices that belong to same segment
        segments = dict()
        for atom in model.atom:
            key = atom.segi
            if key in segments:
                segments[key].append(atom.index)
            else:
                segments[key] = [atom.index]
        network = segments.values()

    else:
        raise CmdException(f"Unknown unwrap mode: {mode!r}")

    # Traverse network and make molecules whole again
    states = _self.count_states(selection)
    for state in range(1, states + 1):
        box = _self.get_symmetry(selection, state)[0:3]
        xyz = _self.get_coords(selection, state)
        for nt in network:
            if len(nt) == 1:
                continue  # Skip if only one atom
            # Apply PBC according to reference (first) atom
            r = xyz[nt]
            dr = r - r[0]
            r -= box * np.rint(dr / box)
            # Check if COM is out of box put it back
            com = np.mean(r, axis=0)
            xyz[nt] = r - box * np.floor(com / box)
        _self.load_coords(xyz, selection, state)


@cmd.extend
def align2eigen(selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with it's eigen vector. 
    USAGE
        align2eigen [ selection [, state ]]
    ARGUMENTS
        selection : str, optional
            Atoms selection. 
        state : int, default = 0
            Which state to consider.
    """
    xyz = _self.get_coords(selection, state)
    # Center coordinates
    mean = np.mean(xyz, axis=0)
    centered = xyz - mean
    # Compute eigenvectors and eigenvalues
    cov = np.cov(centered.T)
    eval, evec = np.linalg.eig(cov)
    # Ensure a proper right-handed coordinate system
    if np.linalg.det(evec) < 0:
        evec[:, -1] *= -1
    # Transform the points to the new basis
    aligned = np.matmul(centered, evec) + mean
    _self.load_coords(aligned, selection, state)


@cmd.extend
def align2points(selection="all", pk1="(pk1)", pk2="(pk2)", pk3="(pk3)", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with three selected point i.e. angle formed by selections 
        `(pk1)`, `(pk2)`, and `(pk3)` which can be set using the `PkAt` mouse action
        (typically, :kbd:`Ctrl` + :kbd:`middle-click`).

        New vector base is defined as::

              y (pk3)
              |
              | 
            (pk2) -- x (pk1)
             /
            z 

    USAGE
        align2points [ selection [, pk1, [, pk2 [, pk3, [, state ]]]]]
    ARGUMENTS
        selection : str, optional
            Atoms selection. 
        pk1, pk2, pk3 : str, default = '(pk1)', '(pk2)', '(pk2)'
            Points that will be used to define new orthonormal base.
        state : int, default = 0
            Which state to consider.
    """
    angle = _self.get_angle(pk1, pk2, pk3)
    if not angle:
        raise ValueError("Selection must be non-colinear points")
    v1 = _self.get_coords(pk1, state)
    v2 = _self.get_coords(pk2, state) # Ref
    v3 = _self.get_coords(pk3, state)
    # Define a new base
    vx = v1 - v2
    vz = np.cross(v3 - v2, vx)
    vy = np.cross(vz, vx)
    # Normalize base
    vx /= np.linalg.norm(vx)
    vy /= np.linalg.norm(vy)
    vz /= np.linalg.norm(vz)
    # Stack to matrix
    base = np.array([vx, vy, vz]).reshape(3, 3).T
    # Transform to new base
    xyz = _self.get_coords(selection, state) - v2
    xyz = np.matmul(xyz, base)
    _self.load_coords(xyz + v2, selection, state=state)
 

@cmd.extend
def split(operator, selection="all", name="obj", *, _self=cmd):
    """
    DESCRIPTION
        Create a single object for each entity in selection, defined by operator
        (e.g. bymolecule, bysegment, ...). Returns the number of created objects.
    USAGE
        split operator [, selection [, name ]]
    ARGUMENTS
        operator : str
            Operator to split selection by (e.g. bymolecule).
        selection : str, optional
            Atoms selection.
        name : str, default = 'obj'
            Name of generated object(s).
    RETURNS
        : int
            Number of created objects.
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    """
    _self.disable(" ".join(_self.get_object_list(selection)))
    tmp = _self.get_unused_name("_")
    _self.create(tmp, selection)

    r = 0
    while _self.count_atoms(tmp) > 0:
        obj = _self.get_unused_name(name)
        _self.extract(obj, operator + " first model " + tmp)
        r += 1

    _self.delete(tmp)
    return r


@cmd.extend
def remove_alt(selection="all", keep="first", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Remove alternative location atoms.
    USAGE
        remove_alt [ selection [, keep ]]
    ARGUMENTS
        selection : str
            Atom selection.
        keep : str | int, default = 'first'
            AltLoc to keep, or `first` to keep the first observed AltLoc
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology 
    """
    if keep == "first":
        alts = {}
        expr = "(alt, q) = callback((model, segi, chain, resi, resn, name), alt)"

        def callback(namekey, alt):
            return ("#", 0.) if alts.setdefault(namekey, alt) != alt else ("", 1.)

        tmpsele = _self.get_unused_name("_sele")
        _self.select(tmpsele, selection)
        try:
            _self.alter(tmpsele, expr, space={"callback": callback})
            _self.remove(f"{tmpsele} & not alt ''", quiet=quiet)
        finally:
            _self.delete(tmpsele)

    else:
        if len(keep) != 1:
            raise CmdException(
                f"keep must be 'first' or a single letter, got {keep!r}")

        _self.remove("(%s) and not alt +%s" % (selection, keep), quiet=int(quiet))
        _self.alter(selection, "(alt,q)=("",1.0)")
        _self.sort()


@cmd.extend
def copy_identifiers(target, source, identifiers="segi chain resi",
                     match="align", *, _self=cmd):
    """
    DESCRIPTION
        Transfers identifiers e.g. segi, chain, and resi from one selection to another.
        This works by mapping old to new identifiers and alters also not aligned
        atoms (works if any other atom from the same residue got aligned).
    USAGE
        copy_identifiers target, source [, identifiers [, match ]]
    ARGUMENTS
        target : str
            Target selection. 
        source : str
            Source selection.
        identifiers : str, default = 'segi chain resi'
            Identifiers to copy. Separator is optional.
        match : str, default = 'align'
            Method how to match atoms.
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    """
    from .fitting import MatchMaker

    with MatchMaker(target, source, match, _self=_self) as mm:
        key = "(" + ",".join(identifiers.split()) + ",)"
        tkeys, skeys = [], []
        _self.iterate(mm.mobile, "tkeys.append(%s)" % (key), space=locals())
        _self.iterate(mm.target, "skeys.append(%s)" % (key), space=locals())
        t2s = dict(zip(tkeys, skeys))
        _self.alter(target, "%s = t2s.get(%s, %s)" %
                    (key, key, key), space=locals())


# Autocomplete
cmd.auto_arg[0].update({
    "renumber": cmd.auto_arg[0]["zoom"],
    "wrap": cmd.auto_arg[0]["zoom"],
    "unwrap": cmd.auto_arg[0]["zoom"],
    "align2eigen": cmd.auto_arg[0]["zoom"],
    "align2points": cmd.auto_arg[0]["zoom"],
    "split": cmd.auto_arg[0]["zoom"],
    "remove_alt": cmd.auto_arg[0]["zoom"],
    "copy_identifiers": cmd.auto_arg[0]["zoom"]
})
cmd.auto_arg[1].update({
    "unwrap": [cmd.Shortcut(_UNWRAP_MODES), "method", ""],
    "copy_identifiers": cmd.auto_arg[0]["zoom"]
})
