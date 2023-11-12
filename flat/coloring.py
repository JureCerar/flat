from pymol import cmd, menu

# Custom color cycle (only 9 colors)
_color_cycle = [
    "0x3ec692", # green
    "0x4c94fc", # blue
    "0xf5a566", # orange
    "0xe05658", # red
    "0xfbda42", # yellow 
    "0xd870ba", # pink
    "0x5c5d5c", # dark
    "0x3fb8c9", # cyan
    "0xaa6fda", # purple
]

@cmd.extend
def cbc(selection="(all)"):
    """
    DESCRIPTION
        Color all chains a different color.
    USAGE
        cbc [ selection ]
    """
    for i, chain in enumerate(cmd.get_chains(selection)):
        if len(chain.split()) != 1:
            chain = f"'{chain}'"
        color = _color_cycle[i % len(_color_cycle)]
        cmd.color(color, f"(({selection}) & chain {chain} & e. C)")
    return


@cmd.extend
def cbe(selection="(all)", *, _self=cmd, **kwargs):
    """
    DESCRIPTION
        Color all elements a different color.
    USAGE
        cbe [ selection [, <element>=<color> ]]
    EXAMPLE
        cbe (all), C=green, oxygen=red
    """
    for element, color in kwargs.items():
        cmd.color(color, f"(({selection}) & e. {element})")
    return


@cmd.extend
def cbattr(selection="(all)", *args, _self=cmd):
    """
    DESCRIPTION
        Color by residue by attribute.
    USAGE
        cbattr [ selection [, attribute ]]
    ATTRIBUTES
        negative, positive, polar, nonpolar, special, element
    """

    if "negative" in args or "all" in args:
        resn = ["ASP", "GLU"]
        _self.color("red", f"({selection}) & resn {'+'.join(resn)}")

    if "positive" in args or "all" in args:
        resn = ["ARG", "HIS", "LYS"]
        _self.color("marine", f"({selection}) & resn {'+'.join(resn)}")

    if "polar" in args or "all" in args:
        resn = ["ASN", "GLN", "SER", "THR"]
        _self.color("tv_green", f"({selection}) & resn {'+'.join(resn)}")

    if "nonpolar" in args or "all" in args:
        resn = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
        _self.color("yellow", f"({selection}) & resn {'+'.join(resn)}")

    if "special" in args or "all" in args:
        resn = ["CYS", "GLY", "PRO"]
        _self.color("gray", f"({selection}) & resn {'+'.join(resn)}")

    if "element" in args or "all" in args:
        _self.color("white", f"({selection}) & e. H")
        _self.color("blue", f"({selection}) & e. N")
        _self.color("red", f"({selection}) & e. O")
        _self.color("yellow", f"({selection}) & e. S")

    return


class _Palette:
    """ Dummy class """
    pass


class _Colors():
    """ Class for storing color palettes. PyMOL API only. """

    def __init__(self, *args):
        self.colors = list(args)
        self.index = 0
        return

    def __repr__(self):
        return " ".join(self.colors)

    def replace(self, *args):
        return " ".join(self.colors)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.colors[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return result

    def list(self):
        return self.colors


# Color palettes. PyMOL API only
palette = _Palette()

# Custom colors
cmd.set_color("flat.green", [62, 198, 146])
cmd.set_color("flat.blue", [76, 148, 252])
cmd.set_color("flat.orange", [245, 165, 102])
cmd.set_color("flat.red", [224, 86, 88])
cmd.set_color("flat.yellow", [251, 218, 66])
cmd.set_color("flat.pink", [216, 112, 186])
cmd.set_color("flat.dark", [92, 93, 92])
cmd.set_color("flat.cyan", [63, 184, 201])
cmd.set_color("flat.purple", [170, 111, 218])
menu.all_colors_list.append(
    ("flat", [
        ("275", "flat.green"),
        ("259", "flat.blue"),
        ("963", "flat.orange"),
        ("833", "flat.red"),
        ("982", "flat.yellow"),
        ("847", "flat.pink"),
        ("333", "flat.dark"),
        ("277", "flat.cyan"),
        ("648", "flat.purple"),
    ]),
)

# Novartis colors
palette.novartis = _Colors("novartis", "carmine", "sienna", "apricot")
palette.novartis_rev = _Colors("apricot", "sienna", "carmine", "novartis")
cmd.set_color("novartis", [4, 96, 169])
cmd.set_color("carmine", [141, 21, 27])
cmd.set_color("sienna", [231, 74, 33])
cmd.set_color("apricot", [236, 154, 30])
menu.all_colors_list.append(
    ('novartis', [
        ('036', 'novartis'),
        ('501', 'carmine'),
        ('921', 'sienna'),
        ('961', 'apricot'),
    ]),
)

# AlphaFold colors
palette.alphafold = _Colors("high", "medium", "low", "crit")
palette.alphafold_rev = _Colors("crit", "low", "medium", "high")
cmd.set_color("high", [0, 83, 214])
cmd.set_color("medium", [101, 203, 243])
cmd.set_color("low", [255, 219, 19])
cmd.set_color("crit", [255, 125, 69])
menu.all_colors_list.append(
    ('alphafold', [
        ('038', 'high'),
        ('379', 'medium'),
        ('980', 'low'),
        ('942', 'crit'),
    ]),
)
