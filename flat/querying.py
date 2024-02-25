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

from pymol import cmd, CmdException, cgo


@cmd.extend
def count(selection="all", *, _self=cmd):
    """
    DESCRIPTION
        Count number of atoms, residues, and mass in selection.
    USAGE
        count [ selection ]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
    """
    # Get properties
    model = _self.get_model(f"({selection})")
    if not model.atom:
        print("Empty selection!")
        return
    mass = model.get_mass()
    print("Residue(s): %s" % len(model.get_residues()))
    print("Atom(s): %s" % model.nAtom)
    if mass < 1e3:
        print("Mass: %.3f Da" % mass)
    elif mass < 1e6:
        print("Mass: %.3f kDa" % (mass/1e3))
    else:
        print("Mass: %.3f MDa" % (mass/1e6))
    return


@cmd.extend
def length(selection="all", num=10, *, _self=cmd):
    """
    DESCRIPTION
        Return the chain, residue, and atoms lengths in selection.
    USAGE
        length [ selection [, num ]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        num = str: Number of list items to display. {default: 10}
    """
    num = int(num)
    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple object in selection.")
    model = _self.get_model(selection)

    def truncate(lst, n):
        """ Return truncated list for writing """
        if n >= len(lst):
            string = str(lst)
        elif n < 1:
            string = str(lst)
        elif n == 1:
            string = str(lst[:1])[:-1]
            string += ", ... ]"
        else:
            lower = int(n // 2)
            upper = len(lst) + lower - n
            string = str(lst[:lower])[:-1]
            string += ", ... "
            string += str(lst[upper:])[1:]
        return string
    
    # Chains
    chain = _self.get_chains(selection)
    print("Chains:", len(chain), truncate(chain, num))
    # Residues
    residue = [int(model.atom[i].resi) for i in range(model.nAtom)]
    residue = list(dict.fromkeys(residue))
    print("Residues:", len(residue), truncate(residue, num))
    # Atoms
    atoms = [model.atom[i].index for i in range(model.nAtom)]
    atoms = list(dict.fromkeys(atoms))
    print("Atoms:", len(atoms), truncate(atoms, num))
    return


@cmd.extend
def get_sasa(selection, state=0, dot_density=4, *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Get solvent accesible surface area for selection.
    USAGE
        get_sasa selection [, state [, dot_density ]]
    """
    state, dot_density, quiet = int(state), int(dot_density), int(quiet)
    if state < 1:
        state = _self.get_state()
    tmp = _self.get_unused_name("_tmp")
    _self.create(tmp, selection, state, 1, zoom=0, quiet=1)
    _self.set("dot_solvent", 1, tmp)
    if dot_density > -1:
        _self.set("dot_density", dot_density, tmp)
    sasa = _self.get_area(tmp, quiet=quiet)
    _self.delete(tmp)
    return sasa


@cmd.extend
def iterate_to_list(selection, expression, *, space=None, _self=cmd):
    """
    DESCRIPTION
        Capture "iterate" results in a list.
    USAGE
        iterate_to_list selection, expression
    """
    outlist = []
    _self.iterate(
        selection,
        "outlist.append(({}))".format(expression),
        space=dict(space or (), outlist=outlist),
    )
    return outlist


@cmd.extend
def ext_coef(selection="all", state=-1, *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Calculate extinction coefficient and absorbance for native and folded protein.
    USAGE
        ext_coef [ selection [, state ]]
    LITERATURE
        C.N. Pace, et. al, Protein Sci., 1995, doi:10.1002/pro.5560041120
        H. Edelhoch, Biochemistry, 1967, doi:10.1021/bi00859a010
        C.G. Gill & P.H. von Hippel, Anal. Biochem., 1989, doi:10.1016/0003-2697(89)90602-7
    """
    state, quiet = int(state), int(quiet)
    nW = _self.count_atoms(f"({selection}) & resn TRP & guide", state=state)
    nY = _self.count_atoms(f"({selection}) & resn TYR & guide", state=state)
    nS = _self.count_atoms(f"({selection}) & resn CYS+CYS2 & guide", state=state)
    nSS = _self.count_atoms(f"({selection}) & (CYS+CYS2/SG and bound_to CYS+CYS2/SG)", state=state)
    # If not higher structure data is avaliable (i.e. calculating from fasta)
    nSS = nS if nSS == 0 else nSS
    # Calculate extinction coefficeint and absorbance at 1 mg/ml for native protein
    mass = _self.get_model(selection).get_mass()
    coef = nW * 5500 + nY * 1490 + nSS * 62.5
    abs = coef / mass
    coef0 = nW * 5500 + nY * 1490
    abs0 = coef0 / mass
    if not quiet:
        print(
            "Extinction coefficients are in units of M^-1 cm^-1, at 280 nm measured in water.",
            "",
            "Ext. coefficient: {:.0f}".format(coef),
            "Abs 0.1% (=1 mg/mL): {:.3f}, for native protein.".format(abs),
            "",
            "Ext. coefficient: {:.0f}".format(coef0),
            "Abs 0.1% (=1 mg/mL): {:.3f}, for denatured protein.".format(abs0),
            sep="\n"
        )
    return ((coef, abs), (coef0, abs0))


@cmd.extend
def get_seq(selection, chainbreak="/", unknown="X", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Gets the one-letter sequence, including residues without coordinates.
    USAGE
        get_seq selection [, chainbreak [, unknown ]]
    """
    # See: https://github.com/speleo3/pymol-psico/blob/master/psico/modelling.py#L427
    from pymol.exporting import _resn_to_aa as one_letter

    seq_list = []
    _self.iterate(
        f"({selection}) & polymer",
        "seq_list.append((resn, resv))",
        space=locals(),
    )

    def seqbuilder():
        prev_resv = None
        for resn, resv in seq_list:
            if resv != prev_resv:
                if prev_resv is not None and resv != prev_resv + 1:
                    yield chainbreak
                if resn in one_letter:
                    yield one_letter[resn]
                else:
                    print('Warning: unknown residue "%s"' % (resn))
                    yield unknown
                prev_resv = resv
    
    seq = seqbuilder()

    if not quiet:
        print("get_seq:", "".join(seq))

    return ''.join(seq)


@cmd.extend
def get_dipole(selection="all", state=0, var="formal_charge", vis=1, quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Get electric dipole momentum for atoms in selection. For VAR use
        either 'partial_charge' or 'formal_charge'. 
    USAGE
        get_dipole [ selection [, state [, var [, vis ]]]]
    ARGUMENTS
        selection = str: Atom selection {default: all}
        state = int: Object state (0 for current state). {default: 0}
        var = str: Property used for calculation. {default: "formal_charge"}
        vis = int: Visualize output. {default: 1}
    """
    import numpy as np
    state, vis, quiet = int(state), int(vis), int(quiet)

    bfact = []
    _self.iterate_state(state, selection, f"bfact.append({var})", space=locals())
    com = np.array(_self.centerofmass(selection, state))
    xyz = np.array(_self.get_coords(selection, state))
    bfact = np.array(bfact, dtype=float)

    if not any(bfact) or len(bfact) == 0:
       raise CmdException(f"Property '{var}' not assigned?")

    dipole = np.zeros(3)
    for i in range(len(xyz)):
        dipole += bfact[i] * (xyz[i] - com)

    dipole /= 0.2081943  # From [e*nm] to [D]

    if not quiet:
        print(
            " Util: Dipole_moment =",
            f"{np.linalg.norm(dipole):.3f} D",
            np.array2string(dipole, precision=2)
        )

    if vis:
        color1 = _self.get_color_tuple("red")
        color2 = _self.get_color_tuple("blue")

        # Scale arrow shape
        radius = 0.5
        hlength = radius * 3.0
        hradius = hlength * 0.6

        vmin, vmax = _self.get_extent(selection, state)
        norm = np.linalg.norm(np.array(vmax) - np.array(vmin))
        v = dipole / np.linalg.norm(dipole) * norm / 2

        xyz1 = com - v
        xyz2 = com + v
        xyz3 = xyz2 + v / np.linalg.norm(v) * hlength

        name = _self.get_unused_name("dipole")
        obj = [
            cgo.CYLINDER, *xyz1, *xyz2, radius, *color1, *color2,
            cgo.CONE, *xyz2, *xyz3, hradius, 0.0, *color2, *color2, 1.0, 0.0,
        ]
        _self.load_cgo(obj, name, state=state, zoom=1)

    return dipole


@cmd.extend
def get_longest_distance(selection="(all)", vis=1, *, _self=cmd):
    """
    DESCRIPTION
        Get longest distance in a selection.
    USAGE
        get_longest_distance [ selection [, vis ]]
    """
    import numpy as np
    from scipy import spatial

    xyz = np.array(_self.get_coords(selection, 1))

    # Find 2 points with largest distance (Convex Hull)
    xyz_hull = xyz[spatial.ConvexHull(xyz).vertices]
    dist_mat = spatial.distance_matrix(xyz_hull, xyz_hull)
    i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
    longest = xyz_hull[i] - xyz_hull[j]

    if vis:
        point_a = _self.get_unused_name("_point_A")
        point_b = _self.get_unused_name("_point_B")
        _self.pseudoatom(point_a, pos=xyz_hull[i].tolist())
        _self.pseudoatom(point_b, pos=xyz_hull[j].tolist())
        _self.distance(
            _self.get_unused_name("distance"),
            point_a,
            point_b,
            mode=0,
        )
        _self.delete(point_a)
        _self.delete(point_b)

    return longest


@cmd.extend
def get_raw_distances(names='', state=1, selection='all', quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Get the list of pair items from distance objects. Each list item is a
        tuple of (index1, index2, distance).
        Based on a script from Takanori Nakane, posted on pymol-users mailing list:
        http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html
    ARGUMENTS
        names = string: names of distance objects (no wildcards!) {default: all measurement objects}
        state = integer: object state {default: 1}
        selection = string: atom selection {default: all}
    SEE ALSO
        select_distances, cmd.find_pairs, cmd.get_raw_alignment
    """
    from chempy import cpv

    state, quiet = int(state), int(quiet)
    if state < 1:
        state = _self.get_state()

    valid_names = _self.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                print(' Error: no such distance object: ' + name)
                raise CmdException

    raw_objects = _self.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    _self.iterate_state(
        state,
        selection,
        'xyz2idx[x,y,z] = (model,index)',
        space=locals()
    )

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append(
                    (
                        xyz2idx[xyz1],
                        xyz2idx[xyz2],
                        cpv.distance(xyz1, xyz2)
                    )
                )
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))
    return r


@cmd.extend
def get_contacts(selection1, selection2, name="contacts", cutoff = [3.6, 4.0, 4.8], state=0, *, _self=cmd):
    """
    DESCRIPTION
        Returns the average number of each contacts between two selections (the good, the bad,
        and the ugly): hydrogen bonds, salt bridges, pi-cation, and pi-pi interactions.
    USAGE
        get_contacts selection1, selection2 [ name [, cutoff [, state ]]]
    ARGUMENTS
        selection1, selection2 = str: Atom selection
        name = str: Name of output group {default: "contacts"}
        cutoff = list: Cutoff distance for each type of interaction {default: [3.6, 4.0, 4.8]}
        state = int: Object state {default: 0 (all states)}
    NOTES
        pi-pi interactions for identical selections not working properly.
    PYTHON API
        cmd.show_contacts(...) -> list
    SEE ALSO
        get_raw_distances
    """
    state = int(state)

    # Check for empty selections
    if not _self.get_object_list(selection1):
        raise CmdException(f"Empty selection: {selection1}")
    if not _self.get_object_list(selection2):
        raise CmdException(f"Empty selection: {selection2}")

    # Generate new empty group
    name = _self.get_legal_name(name)
    _self.delete(name)
    _self.group(name)

    result = [0, 0, 0, 0]

    def get_aromatic(selection, name):
        """ Generate aromatic centers as pseudo-atoms """
        residues = []
        _self.iterate(
            f"({selection}) & (bca. resn PHE+TYR+TRP)",
            "residues.append((chain,resn,resi))",
            space=locals()
        )
        for chain, resn, resi in residues:
            selection = f"///{chain}/{resn}`{resi}/"
            if resn == "PHE":
                selection += "CG+CZ"
            elif resn == "TYR":
                selection += "CG+CZ"
            elif resn == "TRP":
                selection += "CE*+CZ*"
            _self.pseudoatom(name, selection)
        return


    try:
        # Get aromatic pseudo-atoms for each selections
        get_aromatic(selection1, "aromatic1")
        get_aromatic(selection2, "aromatic2")

        # ----------------------
        # hydrogen bonds
        _self.distance(
            "hbond",
            f"{selection1} & (donor|acceptor)",
            f"{selection2} & (donor|acceptor)",
            cutoff=cutoff[0],
            mode=2,
            state=state,
        )
        result[0] = len(get_raw_distances("hbond", state))
        _self.set("dash_color", "yellow", "hbond")
        _self.group(name, "hbond")

        # ----------------------
        # salt bridges
        _self.distance(
            "salt_bridge",
            f"{selection1} & (fc. < 0 | fc. > 0)",
            f"{selection2} & (fc. < 0 | fc. > 0)",
            cutoff=cutoff[1],
            mode=0,
            state=state,
        )
        result[1] = len(get_raw_distances("salt_bridge", state))
        _self.set("dash_color", "red", "salt_bridge")
        _self.group(name, "salt_bridge")

        # ----------------------
        # pi-cation bonds
        _self.distance(
            f"pi-cation",
            f"{selection1} & fc. > 0",
            "?aromatic2",
            cutoff=cutoff[1],
            mode=0,
            state=state,
        )
        _self.distance(
            f"pi-cation",
            "?aromatic1",
            f"{selection2} & fc. > 0",
            cutoff=cutoff[1],
            mode=0,
            state=state,
        )
        result[2] = len(get_raw_distances("pi-cation", state))
        _self.set("dash_color", "purple", "pi-cation")
        _self.group(name, "pi-cation")

        # ----------------------
        # pi-pi bonds
        _self.distance(
            "pi-pi",
            "?aromatic1",
            "?aromatic2",
            cutoff=cutoff[2],
            mode=0,
            state=state,
        )
        result[3] = len(get_raw_distances("pi-pi", state))
        _self.set("dash_color", "gray", "pi-pi")
        _self.group(name, "pi-pi")

    finally:
        _self.delete("aromatic1")
        _self.delete("aromatic2")

    # Normalize to number of frames
    if not state:
        num_states = _self.count_states()
        result = [r / num_states for r in result]

    return result
