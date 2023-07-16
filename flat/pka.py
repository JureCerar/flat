from pymol import cmd
import tempfile
import os

# See: https://propka.readthedocs.io/en/latest/api.html
# See: https://tuilab.github.io/tui/update/2017/02/21/propka.html


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
    import propka.run as pk

    state = int(state)
    vis = bool(vis)

    # Save temporary PDB file
    temp = os.path.join(tempfile.mkdtemp(), "temp.pdb")
    cmd.save(temp, selection, state)

    # Run propka3
    try:
        protein = pk.single(temp, optargs, None, False)
    except:
        pass

    # Get list of residues
    residues = []
    for group in protein.conformations["AVR"].groups:
        if group.residue_type in ["ASP", "GLU", "HIS", "TYRF", "LYS", "ARG", "N+", "C-"]:
            if group.atom.chain_id == "_":
                group.atom.chain_id = ""
            residues.append(group)

    # Visualize atoms and show label
    if vis:
        cmd.show(
            "lines",
            f"((byres ({selection})) & (sc.|(n. CA|n. N&r. PRO)))"
        )
        cmd.set("label_size", 12.0)
        cmd.set("sphere_scale", 0.25, "(all)")
        for group in residues:
            label = "'pKa={:.2f}'".format(group.pka_value)
            atom = "///{}/{}`{}/{}".format(
                group.atom.chain_id,
                group.atom.res_name,
                group.atom.res_num,
                group.atom.name,
            )
            cmd.label(atom, label)
            cmd.show("spheres", atom)

    # Write pKa file
    if file:
        protein.write_pka(file)

    return protein


def get_propka(selection, state=0):
    """
    DESCRIPTION
      Old function
    USAGE
    ARGUMENTS
    """
    import propka.run as pk

    temp = os.path.join(tempfile.mkdtemp(), f"something.pdb")
    cmd.save(temp, selection, state)

    # Run propka3
    try:
        protein = pk.single(temp, [], None, False)
    except:
        pass

    # Get list of residues
    residues = []
    for group in protein.conformations["AVR"].groups:
        if group.residue_type in ["ASP", "GLU", "HIS", "TYR", "LYS", "ARG", "N+", "C-"]:
            if group.atom.chain_id == "_":
                group.atom.chain_id = ""
            residues.append(group)

    # Sort list
    residues = sorted(
        residues,
        key=lambda x: (x.atom.chain_id, x.residue_type)
    )

    # Visualize atoms and show label
    cmd.show_as(
        "lines",
        f"((byres ({selection})) & (sc.|(n. CA|n. N&r. PRO)))"
    )
    cmd.set("label_size", 12.0)
    cmd.set("sphere_scale", 0.25, "(all)")
    for group in residues:
        label = "'pKa={:.2f}'".format(group.pka_value)
        atom = "///{}/{}`{}/{}".format(
            group.atom.chain_id,
            group.atom.res_name,
            group.atom.res_num,
            group.atom.name,
        )
        cmd.label(atom, label)
        cmd.show("spheres", atom)

    # Print residue pKa values
    print("       Group      pKa  model-pKa")
    for group in residues:
        str_ = f"{group.atom.res_name:>6s}"
        str_ += f"{group.atom.res_num:4d}"
        str_ += f"{group.atom.chain_id:>2s}"
        str_ += f"{group.pka_value:9.2f}"
        str_ += f"{group.model_pka:11.2f}"
        print(str_)
    # Print charge profile
    print("-" * 40)
    print("    pH  unfolded  folded")
    charge = protein.get_charge_profile("AVR", [0, 14, 0.5])
    pi = protein.get_pi()
    for point in charge:
        str_ = f"{point[0]:6.1f}"
        str_ += f"{point[1]:10.2f}"
        str_ += f"{point[2]:8.2f}"
        print(str_)
    # Print isoelectric point
    print(f"The pI is {pi[0]:.2f} (folded) and {pi[1]:.2f} (unfolded)")

    return protein
