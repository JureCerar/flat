from pymol import cmd


@cmd.extend
def sum_formal_charges(selection="(all)", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Get sum of formal charges for selection
    USAGE
        sum_formal_charges [selection [, quiet ]]
    """
    _util_sum_fc = [0]
    _self.iterate(
        selection, "_util_sum_fc[0] += formal_charge", space=locals())
    sum_formal_charges = _util_sum_fc[0]
    if not int(quiet):
        print(f" Util: {sum_formal_charges = }")
    return sum_formal_charges


@cmd.extend
def sum_partial_charges(selection="(all)", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Get sum of partial charges for selection
    USAGE
        sum_partial_charges [selection [, quiet ]]
    """
    _util_sum_pc = [0.0]
    _self.iterate(
        selection, "_util_sum_pc[0] += partial_charge", space=locals())
    sum_partial_charges = _util_sum_pc[0]
    if not int(quiet):
        print(f" Util: {sum_partial_charges = }")
    return sum_partial_charges


@cmd.extend
def set_charges_and_radii(obj_name, quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Assign atom VdW radii, formal, and partial charges using AMBER99 force-field.
    USAGE
        set_charges_and_radii obj_name [, quiet ]
    NOTES
        Function has it's own fair share of problems. Use with caution!
    """
    from chempy.champ import assign

    quiet = int(quiet)

    # Convent Seleno-methionine to methionine
    _self.alter(obj_name+"///MSE/SE", "elem='S';name='SD'", quiet=1)
    _self.alter(obj_name+"///MSE/", "resn='MET'", quiet=1)
    _self.flag("ignore", obj_name, "clear")

    # Remove alternate conformers
    _self.remove(obj_name+" and not alt ''+A")
    _self.alter(obj_name, "alt=''")
    _self.sort(obj_name)
    _self.fix_chemistry(obj_name, obj_name, 1)

    # make sure all atoms are included...
    _self.alter(obj_name, "q=1.0", quiet=1)

    if not quiet:
        print(" Util: Fixing termini and assigning formal charges...")

    assign.missing_c_termini(obj_name, quiet=1, _self=_self)

    while not assign.formal_charges(obj_name, quiet=1, _self=_self):
        print(" WARNING: unrecognized or incomplete residues are being deleted:")
        _self.iterate("(byres ("+obj_name+" and flag 23)) and flag 31",
                    'print("  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/")', quiet=1)
        # Get rid of residues that weren't assigned
        _self.remove("byres ("+obj_name+" and flag 23)")
        assign.missing_c_termini(obj_name, quiet=1, _self=_self)

    if not quiet:
        print(" Util: Assigning Amber 99 charges and radii...")

    _self.h_add(obj_name)
    if not assign.amber99(obj_name, quiet=1, _self=_self):
        print(" WARNING: some unassigned atoms are being deleted:")
        _self.iterate("byres ("+obj_name+" and flag 23)",
                    'print("  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/"+name+"? ["+elem+"]")', quiet=1)
        # Get rid of residues that weren't assigned
        _self.remove(obj_name+" and flag 23")

    # Show the user what the net charges are
    formal = sum_formal_charges(obj_name, quiet=0, _self=_self)
    partial = sum_partial_charges(obj_name, quiet=0, _self=_self)
    if round(formal) != round(partial):
        print(" WARNING: formal and partial charge sums don't match -- there is a problem!")


@cmd.extend
def get_vacuum_esp(selection, mode=2, cutoff=10.0, *, _self=cmd):
    """
    DESCRIPTION
        Calculate vacuum electrostatic potential map.
    USAGE
        get_vacuum_esp selection [, mode [, border [, quiet ]]]
    ARGUMENTS
        selection = str: Atom selection.
        mode = int: Calculation mode: {default: 2}
            0: coulomb, no cutoff
            1: coulomb_neutral map, no cutoff
            2: coulomb_local, cutoff
        cutoff = float: Calculation cutoff value {default: 10.0}
    """
    if (selection.split() != [selection] or
            selection not in cmd.get_names('objects')):
        print(" Error: must provide an object name")
        raise cmd.QuietException
    
    obj_name = selection + "_e_chg"
    map_name = selection + "_e_map"
    pot_name = selection + "_e_pot"
    cmd.disable(selection)
    cmd.delete(obj_name)
    cmd.delete(map_name)
    cmd.delete(pot_name)
    cmd.create(obj_name, "((polymer and ("+selection +
        ") and (not resn A+C+T+G+U)) or ((bymol (polymer and (" +
        selection+"))) and resn NME+NHE+ACE)) and (not hydro)"
    )
    
    # try to just get protein...
    set_charges_and_radii(obj_name, _self=_self)

    ext = cmd.get_extent(obj_name)
    max_length = max(
        abs(ext[0][0] - ext[1][0]),
        abs(ext[0][1] - ext[1][1]),
        abs(ext[0][2] - ext[1][2])
    ) + 2 * cutoff

    # compute an grid with a maximum dimension of 50, with 10 A borders around molecule, and a 1.0 A minimum grid
    sep = max_length/50.0
    if sep < 1.0:
        sep = 1.0

    print(" Util: Calculating electrostatic potential...")
    if mode == 0:  # absolute, no cutoff
        cmd.map_new(map_name, "coulomb", sep, obj_name, cutoff)
    elif mode == 1:  # neutral, no cutoff
        cmd.map_new(map_name, "coulomb_neutral", sep, obj_name, cutoff)
    else:  # local, with cutoff
        cmd.map_new(map_name, "coulomb_local", sep, obj_name, cutoff)

    # Visualize result
    cmd.ramp_new(pot_name, map_name, selection=obj_name, zero=1)
    cmd.hide("everything", obj_name)
    cmd.show("surface", obj_name)
    cmd.set("surface_color", pot_name, obj_name)
    cmd.set("surface_ramp_above_mode", 1, obj_name)