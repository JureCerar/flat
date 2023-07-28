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
    
# @cmd.extend
# def save_fasta( filename, selection='(all)', invert=0, gapped=1, quiet=1 ):
#     '''
#     Synopsis
#         saveFasta filename [, selection [, invert [, gapped]]]
#     Desciption:
#         Save a html file with colored (by C-alpha atoms) fasta sequence.
#     Source:
#         save_colored_fasta(): https://github.com/speleo3/pymol-psico/blob/master/psico/fasta.py
#     Arguments:
#         filename    ...
#         selection   Object selection (default=all).
#         invert      Invert foreground and background colors (default=0)
#         gapped      ...
#     Example
#         saveFasta fasta.html
#     '''
#     from pymol.exporting import _resn_to_aa as one_letter
#     from pymol import Scratch_Storage
#     one_letter.update( {'CYS2':'C','ASPP':'D','GLUP':'E','HSE':'H','HSD':'H','HSP':'H','LSN':'K'} )
#     gapped, invert= int(gapped), int(invert)
#     selection = '(%s) and guide' % (selection)
#     html = []
#     stored = Scratch_Storage()
#     def callback(resv, resn, color):
#         if stored.resv is None:
#             stored.resv = resv - (resv % 70)
#         if gapped:
#             while stored.resv+1 < resv:
#                 callback(stored.resv+1, '-', 25)
#         stored.resv += 1
#         if stored.resv % 70 == 1:
#             html.append(('</span>\n<br>%4d <span>' % (resv)).replace(' ', '&nbsp;'))
#             stored.color = None
#         rgb = cmd.get_color_tuple(color)
#         rgb = [ int(x*255) for x in rgb ]
#         color = '#%02x%02x%02x' % (rgb[0],rgb[1],rgb[2])
#         if not invert:
#             # https://stackoverflow.com/questions/11867545    
#             textColor =  ['white','black'][ round((rgb[0]*299+rgb[1]*587+rgb[2]*114)/1000) > 125 ]
#             backgroundColor = color
#         else:
#             textColor = color
#             backgroundColor='white'
#         aa = one_letter.get(resn,'-')
#         if color != stored.color:
#             html.append('</span><span style="color:'+textColor+';background-color:'+backgroundColor+'">')
#             stored.color = color
#         html.append(aa)
#         if stored.resv % 10 == 0:
#             html.append("</span> <span>")
#     for obj in cmd.get_object_list(selection):
#         for chain in cmd.get_chains('model %s and (%s)' % (obj, selection)):
#             sele = 'model %s and chain "%s" and (%s)' % (obj, chain, selection)
#             html.append('\n<br>&gt;%s_%s<span>' % (obj, chain))
#             stored.resv = None if gapped else 0
#             stored.color = None
#             cmd.iterate(sele, 'callback(resv, resn, color)', space=locals())
#             html.append('</span>')
#     handle = open( filename, 'w')
#     print( '<html><body style="font-family:consolas">', file=handle )
#     print( ''.join(html), file=handle )
#     print( '</body></html>', file=handle )
#     handle.close()
#     return 
