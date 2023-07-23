from pymol import cmd
import collections
import re


@cmd.extend
def save_bfact(file, selection="(all)", var="b", *, _self=cmd):
    """
    DESCRIPTION
        Save property from selection to CSV file.
    USAGE
        saveBfact file [, selection [, var ]]
    """
    bfact = []
    _self.iterate(
        selection,
        "bfact.append([model,segi,chain,resn,resi,name,{}])".format(var),
        space=locals(),
    )
    with open(file, "w") as handle:
        for b in bfact:
            print(*b, sep=",", file=handle)
    return


@cmd.extend
def load_bfact(file, selection="(all)", var="b", *, _self=cmd):
    """
    DESCRIPTION
        Load property to selection from CSV file.
    USAGE
        loadBfact file [, selection [, var ]]
    """
    vis = int(vis)
    bfact = collections.defaultdict()
    for line in open(file, "r"):
        col = line.split(",")
        key, val = col[:-1], col[-1]
        bfact[key] = val
    _self.alter(
        selection,
        "{} = bfact[model,segi,chain,resn,resi,name]".format(var),
        space=locals()
    )
    return


def _get_keyword(string):
    """ Scan a string for keywords i.e. '[ keyword ]' """
    k = re.search(r'\[(.*?)\]', string)
    return k.string.split()[1] if k else None


@cmd.extend
def load_topol(file, selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Load GROMACS topology (.top) file to selection.
    USAGE
        load_topol file [, selection ]
    """
    atom_list = collections.defaultdict(lambda: [0, 0])
    current = None

    for line in open(file, "r"):
        line = line.split(";")[0].strip()
        if not line:
            continue

        key = _get_keyword(line)
        if key:
            current = key
            continue

        if current == "atoms":
            # nr, type, resnr, residue, atom, cgnr, pcharge, mass
            col = line.split()
            try:
                index, pcharge, mass = int(col[0]), float(col[6]), float(col[7])
                atom_list[index] = [pcharge, mass]
            except:
                print(f"Cannot parse line: '{line}'")

    _self.alter(
        selection,
        "partial_charge, mass = atom_list[index]",
        space=locals()
    )
    return


@cmd.extend
def save_ndx(file="index.ndx", quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Save all selections as GROMACS index (.ndx) file.
    USAGE
        save_ndx [ file [, quiet ]]
    """
    quiet = int(quiet)

    selections = cmd.get_names("public_selections")
    groups = dict()

    for name in selections:
        groups[name] = list()
        cmd.iterate(name, "group.append(index)", space={"group": groups[name]})

    with open(file, "w") as handle:
        nlines = 15  # Number of indeces per line
        for name, group in groups.items():
            if not quiet:
                print(f"Saving group: '{name}' len(group atoms)")
            print(f"[ {name} ]", file=handle)
            for i in range(0, len(group), nlines):
                print(" ".join(f"{x:5}" for x in group[i:i+nlines]), file=handle)

    return


@cmd.extend
def load_ndx(file, quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Read a GROMACS index (.ndx) file as a selection.
    USAGE
        read_ndx file [, quiet ]
    """
    quiet = int(quiet)

    groups = dict()

    current = None
    for line in open(file, "r"):
        line = line.split(";")[0]
        key = _get_keyword(line)
        if key:
           current = key
           groups[current] = []
        else:
           groups[current] += [int(i) for i in line.split()]

    for name, group in groups.items():
        # Pymol does not like '&' and '|' characters and reserved keywords
        name = name.replace("&", "_and_")
        name = name.replace("|", "_or_")
        if name.lower() in ["sidechain", "backbone"]:
            name = name + "_"

        if len(group) == 0:
            print(" Empty group: '%s'" % name)
            continue

        if not quiet:
            print("Loading group: '{}' ({} atoms)".format(name, len(group)))
        
        # Pymol does not like long selection strings so we select the group iteratively
        every = 15
        cmd.select(f"({name})", f"index {group[0]}")
        for i in range(1, len(group), every):
            buffer = "index " + "+".join(str(x) for x in group[i:i+every])
            cmd.select(name, f"({name}) | {buffer}", quiet=1)

    return







