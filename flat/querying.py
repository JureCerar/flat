from pymol import cmd, CmdException, cgo


@cmd.extend
def count(selection="all"):
    """
    DESCRIPTION
        Count number of atoms, residues, and mass in selection.
    USAGE
        count [ selection ]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
    """
    # Get properties
    model = cmd.get_model(f"({selection})")
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
def length(selection="all", num=10):
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
    if len(cmd.get_object_list(selection)) > 1:
        raise CmdException("Multiple object in selection.")
    model = cmd.get_model(selection)

    def truncate(lst, n):
        """ Return truncated list for writing. """
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
    chain = cmd.get_chains(selection)
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
        state = cmd.get_state()
    tmp = cmd.get_unused_name("_tmp")
    cmd.create(tmp, selection, state, 1, zoom=0, quiet=1)
    cmd.set("dot_solvent", 1, tmp)
    if dot_density > -1:
        cmd.set("dot_density", dot_density, tmp)
    sasa = cmd.get_area(tmp, quiet=quiet)
    cmd.delete(tmp)
    return sasa


@cmd.extend
def iterate_to_list(selection, expression, *, space=None):
    """
    DESCRIPTION
        Capture "iterate" results in a list.
    USAGE
        iterate_to_list selection, expression
    """
    outlist = []
    cmd.iterate(
        selection,
        "outlist.append(({}))".format(expression),
        space=dict(space or (), outlist=outlist),
    )
    return outlist


@cmd.extend
def ext_coef(selection="all", state=-1, *, quiet=0):
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
    nW = cmd.count_atoms(f"({selection}) & resn TRP & guide", state=state)
    nY = cmd.count_atoms(f"({selection}) & resn TYR & guide", state=state)
    nS = cmd.count_atoms(f"({selection}) & resn CYS+CYS2 & guide", state=state)
    nSS = cmd.count_atoms(
        f"({selection}) & (CYS+CYS2/SG and bound_to CYS+CYS2/SG)",
        state=state
    )
    # If not higher structure data is avaliable (i.e. calculating from fasta)
    nSS = nS if nSS == 0 else nSS
    # Calculate extinction coefficeint and absorbance at 1 mg/ml for native protein
    mass = cmd.get_model(selection).get_mass()
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
def get_seq(selection, chainbreak="/", unknown="X", *, _self=cmd):
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
        from . import one_letter
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

    return ''.join(seqbuilder())


@cmd.extend
def get_dipole(selection="all", state=0, var="partial_charge", vis=1, *, _self=cmd):
    """
    DESCRIPTION
        Get electric dipole momentum for atoms in selection. For VAR use
        either 'partial_charge' or 'formal_charge'. 
    USAGE
        get_dipole [ selection [, state [, var [, vis ]]]]
    ARGUMENTS
        selection = str: Atom selection {default: all}
        state = int: Object state (0 for current state). {default: 0}
        var = str: Propedrty used for calculation. {default: "partial_charge"}
        vis = int: Visualize output. {default: 1}
    """
    import numpy as np
    state, vis = int(state), int(vis)

    bfact = []
    cmd.iterate_state(state, selection, f"bfact.append({var})", space=locals())
    com = np.array(cmd.centerofmass(selection, state))
    xyz = np.array(cmd.get_coords(selection, state))

    if not any(bfact) or len(bfact) == 0:
       raise CmdException(f"Property '{var}' not assigned?")

    dipole = np.zeros(3)
    for i in range(len(xyz)):
        dipole += bfact[i] * (xyz[i] - com)

    dipole /= 0.2081943  # From [e*nm] to [D]

    if vis:
        color1 = cmd.get_color_tuple("red")
        color2 = cmd.get_color_tuple("blue")

        radius = 0.5
        hlength = radius * 3.0
        hradius = hlength * 0.6

        vmin, vmax = cmd.get_extent(selection, state)
        norm = np.linalg.norm(np.array(vmax) - np.array(vmin))
        v = dipole / np.linalg.norm(dipole) * norm / 2

        xyz1 = com - v
        xyz2 = com + v
        xyz3 = xyz2 + v / np.linalg.norm(v) * hlength

        name = cmd.get_unused_name("dipole_moment")
        obj = [cgo.CYLINDER, *xyz1, *xyz2, radius, *color1, *color2]
        obj += [cgo.CONE, *xyz2, *xyz3, hradius,
                0.0, *color2, *color2, 1.0, 0.0]

        cmd.load_cgo(obj, name, state=state, zoom=1)

    return dipole
