from pymol import cmd, CmdException
import collections
import re


@cmd.extend
def load_csv(filename, selection="(all)", var="b", vis=0, *, _self=cmd):
    """
    DESCRIPTION
        Load property from CSV file to selection.
    USAGE
        load_csv filename [, selection [, var [, vis ]]]
    ARGUMENTS
        filename = str: Input file name.
        selection = str: Atom selection. {default: all}
        var = str: Property to save. {default: b}
        vis = int: Visualize output. {default: 1}
    """
    vis = bool(vis)
    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple objects in selection.")

    data = collections.defaultdict(float)
    with open(filename, "r") as handle:
        # From header try to guess CSV format
        col = handle.readline().split(",")
        if len(col) == 5:
            mode = "atom"
        elif len(col) == 4:
            mode = "res"
        else:
            raise CmdException("Cannot determine CSV file format from header")

        # Read the data from the rest of the file
        for line in handle:
            col = line.strip().split(",")
            key, val = tuple(col[:-1]), col[-1]
            data[key] = float(val)

    if mode == "atom":
        # Set property by atom
        _self.alter(
            selection,
            "{} = data[chain, resn, resi, name]".format(var),
            space=locals(),
        )
    elif mode == "res":
        # Set property by residues
        _self.alter(
            selection,
            "{} = data[chain, resn, resi]".format(var),
            space=locals(),
        )

    # Visualize the property
    if vis:
        palette = ["marine", "silver", "red"]
        obj = _self.get_object_list(selection)[0]
        range = (min(data.values()), max(data.values()))
        _self.spectrum(var, " ".join(palette), selection, range[0], range[1])
        _self.ramp_new("ramp", obj, range, palette)

    return


@cmd.extend
def load_topol(filename, selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Load GROMACS topology (.top) file to selection.
    USAGE
        load_topol filename [, selection ]
    """
    atom_list = collections.defaultdict(lambda: [0, 0])
    current = None

    def get_keyword(string):
        """ Scan a string for keywords i.e. '[ keyword ]' """
        k = re.search(r'\[(.*?)\]', string)
        return k.string.split()[1] if k else None

    for line in open(filename, "r"):
        line = line.split(";")[0].strip()
        if not line:
            continue

        key = get_keyword(line)
        if key:
            current = key
            continue

        if current == "atoms":
            # nr, type, resnr, residue, atom, cgnr, pcharge, mass
            col = line.split()
            try:
                index, pcharge, mass = int(
                    col[0]), float(col[6]), float(col[7])
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
def load_ndx(filename, quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Read a GROMACS index (.ndx) file as a selection.
    USAGE
        read_ndx filename [, quiet ]
    """
    quiet = int(quiet)
    groups = dict()

    def get_keyword(string):
        """ Scan a string for keywords i.e. '[ keyword ]' """
        k = re.search(r'\[(.*?)\]', string)
        return k.string.split()[1] if k else None

    current = None
    for line in open(filename, "r"):
        line = line.split(";")[0]
        key = get_keyword(line)
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
        _self.select(f"({name})", f"index {group[0]}")
        for i in range(1, len(group), every):
            buffer = "index " + "+".join(str(x) for x in group[i:i+every])
            _self.select(name, f"({name}) | {buffer}", quiet=1)

    return


@cmd.extend
def save_colored_fasta(file, selection="all", invert=0, *, _self=cmd):
    """
    DESCRIPTION
        Save a HTML file with colored (by C-alpha atoms) fasta sequence.
    USAGE
        save_colored_fasta file [, selection [, invert ]]]
    """
    # See: https://github.com/speleo3/pymol-psico/blob/master/psico/fasta.py
    from pymol import Scratch_Storage
    from pymol.exporting import _resn_to_aa as one_letter

    invert = int(invert)
    selection = f"({selection}) and guide"

    html = []
    stored = Scratch_Storage()

    def callback(resv, resn, color):
        """ Translate to residue to HTML."""
        if stored.resv is None:
            stored.resv = resv - (resv % 70)
        while stored.resv+1 < resv:
            callback(stored.resv+1, None, 222)
        stored.resv += 1
        if stored.resv % 70 == 1:
            html.append("<br>")

        rgb = cmd.get_color_tuple(color)
        rgb = [int(x*255) for x in rgb]
        aa = one_letter.get(resn, "?") if resn else "-"

        if not resn:
            color = "#000000"
            background = ""
        elif invert:
            cs = round((rgb[0] * 299 + rgb[1] * 587 + rgb[2] * 114) / 1000)
            color = "black" if cs > 125 else "white"
            background = "background-color:#{:02x}{:02x}{:02x}".format(*rgb)
        else:
            color = "#{:02x}{:02x}{:02x}".format(*rgb)
            background = ""

        html.append(
            f"<span style=color:{color};{background}>{aa}</span>"
        )
        if stored.resv % 10 == 0:
            html.append("<span> </span>")
        return

    for obj in _self.get_object_list(selection):
        for chain in _self.get_chains(f"model {obj} and ({selection})"):
            sele = f"m. {obj} & c. '{chain}' and ({selection})"
            html.append(f"&gt;{obj}|{chain}")
            stored.resv = 0
            _self.iterate(sele, "callback(resv, resn, color)", space=locals())
            html.append("<br>")

    with open(file, "w") as handle:
        print("<html><body style='font-family:consolas'>", file=handle)
        print("".join(html), file=handle)
        print("</body></html>", file=handle)

    return


try:
    from pymol.importing import loadfunctions
    loadfunctions.setdefault("csv", load_csv)
    loadfunctions.setdefault("top", load_topol)
    loadfunctions.setdefault("ndx", load_ndx)
except ImportError:
    pass
