from pymol import cmd, CmdException
import collections
import sys


@cmd.extend
def renumber(selection, start=1, quiet=False):
    """
    DESCRIPTION
        Renumber sets new residue numbers (resi) for a polymer based on connectivity. 
    ARGUMENTS
        selection = string: atom selection to renumber {default: all}
        start = integer: counting start {default: 1}
    """
    # See: https://pymolwiki.org/index.php/Renumber
    start, quiet = int(start), bool(quiet)
    model = cmd.get_model(selection)
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(10**5)
    cmd.iterate(
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
    cmd.alter(
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
    DESCRIPTON
        Relax the given selection.
    USAGE
        relax selection [, backbone [, neighbors [, cycles [, state ]]]]
    """
    backbone, neighbors,  = int(backbone), int(neighbors)
    cycles, state = int(cycles), int(state)

    obj = cmd.get_object_list(selection)[0]

    cmd.protect()
    cmd.deprotect(selection)
    if not backbone:
        cmd.protect("name CA+C+N+O+OXT")
    if neighbors:
        cmd.deprotect(f"byres ({obj} within 6.0 of ({selection}))")

    cmd.sculpt_activate(obj, state)

    cmd.set("sculpt_vdw_weight", 0.25, obj)  # Low VDW forces
    cmd.set("sculpt_field_mask", 0x1FF, obj)  # Default
    cmd.sculpt_iterate(obj, state, int(cycles * 0.75))

    cmd.set("sculpt_field_mask", 0x01F, obj)  # Local Geometry Only
    cmd.sculpt_iterate(obj, state, int(cycles * 0.25))

    cmd.unset("sculpt_vdw_weight", obj)
    cmd.unset("sculpt_field_mask", obj)
    cmd.sculpt_deactivate(obj)
    cmd.deprotect()

    return


@cmd.extend
def graft(mobile, target, out=None, sculpt=1, cycles=100, *, _self=cmd):
    """
    DESCRIPTION
        Graft mobile structure on target structure. Similarity between aligned
        mobile and target structure must be >80 %.
    USAGE
        graft mobile, target [, out [, minimize ]]
    ARGUMENTS
        mobile = string: atom selection of mobile section.
        target = string: atom selection of target section. 
        out = string: output structure name {default: None}
        minimize = str: Do a short minimization of grafted structure {default: True}
    """
    sculpt, cylces = str(sculpt), str(cycles)

    if len(cmd.get_chains(mobile)) != 1:
        raise CmdException("Mobile selection must contain only one chain")

    if len(cmd.get_chains(target)) != 1:
        raise CmdException("Target selection must contain only one chain")

    _mobile = cmd.get_unused_name("mobile")
    cmd.copy_to(_mobile, mobile, rename="", zoom=0)
    mobile_len = len(cmd.get_model(_mobile).get_residues())

    _target = cmd.get_unused_name("target")
    cmd.copy_to(_target, target, rename="", zoom=0)
    target_len = len(cmd.get_model(_target).get_residues())

    # Clean-up mobile section; chain and segi must match target section.
    model = cmd.get_model(_target)
    chain, segi = model.atom[0].chain, model.atom[0].segi
    cmd.alter(_mobile, f"chain='{chain}'; segi='{segi}'")
    cmd.remove(f"{_mobile} & name OXT")

    # Aligment method "cealign" method works best for this use
    cmd.extra_fit(_mobile, _target, "cealign",
                  mobile_state=-1, target_state=-1)
    align = cmd.get_unused_name("align")
    score = cmd.align(_mobile, _target, quiet=0, object=align)
    if score[0] > 100.0:
        raise Warning("RMSD value is high. Are you sure this is correct aligment?")

    # Renumber target and mobile to correct numbering
    model = cmd.get_model(f"{align} & {_target}")
    first, last = int(model.atom[0].resi), int(model.atom[-1].resi)

    if not out:
        out = cmd.get_unused_name("out")
    renumber(_mobile, start=first)
    cmd.copy_to(out, _mobile, rename="", zoom=0)

    # Split target section into left and right objects.
    # Correctly renumber objects and then combine/join them.
    # Do a energy minimization around newly bonded sections.
    if first > 1:
        target_left = cmd.get_unused_name("target_left")
        cmd.select(target_left, f"{_target} & resi 1-{first-1}")
        renumber(target_left, start=1)
        cmd.copy_to(out, target_left, rename="", zoom=0)
        cmd.bond(
            f"/{out}///{first-1}/C",
            f"/{out}///{first}/N",
        )
        cmd.delete(target_left)
        if sculpt:
            relax(f"{out} & resi {first}-{first-1}", 1, 1, cycles)

    if last < target_len:
        target_right = cmd.get_unused_name("target_right")
        cmd.select(target_right, f"{_target} & resi {last+1}-{target_len}")
        nbond = first+mobile_len
        renumber(target_right, start=nbond)
        cmd.copy_to(out, target_right, rename="", zoom=0)
        cmd.bond(
            f"/{out}///{nbond-1}/C",
            f"/{out}///{nbond}/N",
        )
        cmd.delete(target_right)
        if sculpt:
            relax(f"{out} & resi {nbond-1}-{nbond}", 1, 1, cycles)

    # Clean-up objects (in order they were created)
    cmd.delete(_mobile)
    cmd.delete(_target)
    cmd.delete(align)

    return


@cmd.extend
def mutate(selection, residue, sculpt=0, cycles=100, *, _self=cmd):
    """
    DESCRIPTION
        Mutate a single residue and select best rotamer. Can do some
        scultiong at the end to best rotamer.
    USAGE
        mutate selection, residue [, sculpt [, cycles ]]
    """
    sculpt = int(sculpt)
    
    if len(cmd.get_model(selection).get_residues()) != 1:
        raise CmdException("Multiple residues in selection.")

    try:
        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode(residue.upper())
        cmd.get_wizard().do_select(f"byres ({selection})")

        scores = cmd.get_wizard().bump_scores
        if scores:
            frame = scores.index(min(scores)) + 1
        else:
            frame = 1
        cmd.frame(str(frame))
        cmd.get_wizard().apply()

    finally:
        cmd.set_wizard()

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

    if len(cmd.get_model(selection).get_residues()) != 1:
        raise CmdException("Multiple residues in selection.")
    if len(cmd.get_model(fragment).get_residues()) != 1:
        raise CmdException("Multiple residues in fragment.")

    temp = cmd.get_unused_name("_tmp")
    cmd.create(temp, fragment, source_state=1, target_state=1)

    model = cmd.get_object_list(selection)[0]
    at = cmd.get_model(selection).atom[0]
    segi, chain, resi = at.segi, at.chain, int(at.resi)
    cmd.alter(temp, f"segi,chain,resi='{segi}','{chain}',{resi}")

    cmd.align(
        mobile=f"({temp}) & n. O+C+N",
        target=f"({selection}) & n. O+C+N",
    )

    cmd.remove(selection)
    cmd.fuse(model, temp, mode=3)
    cmd.edit()

    cmd.bond(
        atom1=f"/{model}/{segi}/{chain}/`{resi-1}/C",
        atom2=f"/{model}/{segi}/{chain}/`{resi}/N",
    )
    cmd.bond(
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
    """
    # See: https://github.com/speleo3/pymol-psico/blob/master/psico/modelling.py#L249
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
    cmd.iterate(
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
                cmd.wizard("mutagenesis")
                cmd.do("refresh_wizard")
                cmd.get_wizard().set_mode(resn)
                sele = "/{}/{}/{}/{}`{}".format(*key)
                cmd.get_wizard().do_select(sele)
                scores = cmd.get_wizard().bump_scores
                if scores:
                    frame = scores.index(min(scores)) + 1
                else:
                    frame = 1
                cmd.frame(str(frame))
                cmd.get_wizard().apply()
                sele_list.appennd(sele)
                # if not quiet:
                #     print(f"Mutated: {sele}")

            except Exception as ex:
                print(f"Mutating '{sele}' failed: {ex}")

            finally:
                cmd.set_wizard()

    if sculpt:
        for sele in sele_list:
            relax(sele, 0, 0, cycles)

    return 
