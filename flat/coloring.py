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
        cmd.color(color, f"(({selection}) & {chain} & e. C)")
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
def cbattr(selection="(all)", negative=0, positive=0, polar=0, nonpolar=0, special=0, *, _self=cmd):
    """
    DESCRIPTION
        Color by residue by attribute.
    USAGE
        cbattr [ selection [, attribute ]]
    """
    negative = int(negative)
    positive = int(positive)
    polar = int(polar)
    nonpolar = int(nonpolar)
    special = int(special)

    if not any([negative, positive, polar, nonpolar, special]):
        negative, positive, polar, nonpolar, special = 1, 1, 1, 1, 1

    _res_negative = ["ASP", "GLU"]
    _res_positive = ["ARG", "HIS", "LYS"]
    _res_polar = ["ASN", "GLN", "SER", "THR"]
    _res_nonpolar = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
    _res_special = ["CYS", "GLY", "PRO"]

    if negative:
        cmd.color("red", f"({selection}) & resn {'+'.join(_res_negative)}")
    if positive:
        cmd.color("marine", f"({selection}) & resn {'+'.join(_res_positive)}")
    if polar:
        cmd.color("tv_green", f"({selection}) & resn {'+'.join(_res_polar)}")
    if nonpolar:
        cmd.color("yellow", f"({selection}) & resn {'+'.join(_res_nonpolar)}")
    if special:
        cmd.color("gray", f"({selection}) & resn {'+'.join(_res_special)}")

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
cmd.set_color("my_green", [62, 198, 146])
cmd.set_color("my_blue", [76, 148, 252])
cmd.set_color("my_orange", [245, 165, 102])
cmd.set_color("my_red", [224, 86, 88])
cmd.set_color("my_yellow", [251, 218, 66])
cmd.set_color("my_pink", [216, 112, 186])
cmd.set_color("my_dark", [92, 93, 92])
cmd.set_color("my_cyan", [63, 184, 201])
cmd.set_color("my_purple", [170, 111, 218])
menu.all_colors_list.append(
    ("custom", [
        ("275", "my_green"),
        ("259", "my_blue"),
        ("963", "my_orange"),
        ("833", "my_red"),
        ("982", "my_yellow"),
        ("847", "my_pink"),
        ("333", "my_dark"),
        ("277", "my_cyan"),
        ("648", "my_purple"),
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
palette.alphafold = _Colors("af_high", "af_medium", "af_low", "af_crit")
palette.alphafold_rev = _Colors("af_crit", "af_low", "af_medium", "af_high")
cmd.set_color("af_high", [0, 83, 214])
cmd.set_color("af_medium", [101, 203, 243])
cmd.set_color("af_low", [255, 219, 19])
cmd.set_color("af_crit", [255, 125, 69])
menu.all_colors_list.append(
    ('alphafold', [
        ('038', 'af_high'),
        ('379', 'af_medium'),
        ('980', 'af_low'),
        ('942', 'af_crit'),
    ]),
)
