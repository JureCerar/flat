from pymol import cmd
import tempfile
import os

@cmd.extend
def propka(selection="all", state=0, file=None, vis=1, optargs=[]):
    """
    DESCRIPTION
        Predicts the pKa values of ionizable groups in proteins and
        protein-ligand complexes based on the 3D structure using PROPKA.
    USAGE
        propka [ selection [, state [, file [, vis [, optargs ]]]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        state = int: Object state (0 for all states). {default: 0}
        file = str: Save PROPKA output to file. {default: None}
        vis = bool: Visualize residue pKa values. {default: True}
        optargs = list: Optional arguments to be passed to PROPKA. {default: []}
    """
    # See: https://propka.readthedocs.io/en/latest/api.html
    # See: https://tuilab.github.io/tui/update/2017/02/21/propka.html
    import propka.run as pk

    state, vis = int(state), int(vis)

    temp = os.path.join(tempfile.mkdtemp(), "temp.pdb")
    cmd.save(temp, selection, state)

    try:
        pka = pk.single(temp, optargs, None, False)
    except:
        pass

    if vis:
        cmd.show(
            "lines",
            f"((byres ({selection})) & (sc.|(n. CA|n. N&r. PRO)))"
        )
        cmd.set("label_size", 12.0)
        cmd.set("sphere_scale", 0.25, "(all)")
        for group in pka.conformations["AVR"].groups:
            if group.residue_type in ["ASP", "GLU", "HIS", "TYRF", "LYS", "ARG", "N+", "C-"]:
                chain_id = "" if group.atom.chain_id == "_" else group.atom.chain_id
                label = "'pKa={:.2f}'".format(group.pka_value)
                atom = "///{}/{}`{}/{}".format(
                    chain_id,
                    group.atom.res_name,
                    group.atom.res_num,
                    group.atom.name,
                )
                cmd.label(atom, label)
                cmd.show("spheres", atom)

    if file:
        pka.write_pka(file)

    return pka

@cmd.extend
def anarci(selection="all"):
    return