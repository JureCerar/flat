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
import numpy as np


@cmd.extend
def align_eig(selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with along it's eigen vector. 
    USAGE
        align_eig [ selection [, state ]]
    """
    xyz = _self.get_coords(selection, state)
    cov = np.cov(xyz, rowvar=False)
    eval, evec = np.linalg.eig(cov)
    xyz -= np.mean(xyz, axis=0)
    if np.linalg.det(evec): 
        evec = np.linalg.inv(evec)
    xyz = np.dot(xyz, evec)
    _self.load_coords(xyz, selection, state=state)
    return


@cmd.extend
def align_3p(pk1="pk1", pk2="pk2", pk3="pk3", selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with along three selected points. 
    USAGE
        align_3p pk1, pk2, pk3, [ selection [, state ]]
    """
    v1 = _self.get_coords(pk1, state)
    v2 = _self.get_coords(pk2, state) 
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
    base = np.array([vx, vy, vz]).reshape(3, 3)
    # Transform to new base
    # TODO: Account for center shift
    xyz = _self.get_coords(selection, state).T
    xyz = np.matmul(base.T, xyz)
    _self.load_coords(xyz.T, selection, state=state)
    return


@cmd.extend
def split(operator, selection, prefix="entity", *, _self=cmd):
    """
    DESCRIPTION
        Create a single object for each entity in selection, defined by operator
        (e.g. bymolecule, bysegment, ...). Returns the number of created objects.
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    """
    _self.disable(" ".join(_self.get_object_list(selection)))
    tmp = _self.get_unused_name("_")
    _self.create(tmp, selection)

    r = 0
    while _self.count_atoms(tmp) > 0:
        name = _self.get_unused_name(prefix)
        _self.extract(name, operator + " first model " + tmp)
        r += 1

    _self.delete(tmp)
    return r


@cmd.extend
def remove_alt(selection="all", keep="first", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Remove alternative location atoms.
    USAGE
        remove_alt [selection [, keep]]
    ARGUMENTS
        selection = string: atom selection
        keep = string: AltLoc to keep, or "first" to keep the first observed AltLoc {default: first}
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

    return


@cmd.extend
def copy_identifiers(target, source, identifiers="segi chain resi",
                       match="align", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Transfers identifiers e.g. segi, chain, and resi from one selection to another.
        This works by mapping old to new identifiers and alters also not aligned
        atoms (works if any other atom from the same residue got aligned).
    USAGE
        copy_identifiers target, source [, identifiers [, match ]]
    ARGUMENTS
        target = string: target selection 
        source = string: source selection
        identifiers = string: identifiers to copy {default: "segi,chain,resi"}
        match = string: method how to match atoms {default: "align"}
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
    "split": cmd.auto_arg[0]["zoom"],
    "remove_alt": cmd.auto_arg[0]["zoom"],
    "copy_identifiers": cmd.auto_arg[0]["zoom"]
})
