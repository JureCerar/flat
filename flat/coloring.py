from pymol import cmd, menu

# hould match the list in Color.cpp
# See: https://github.com/schrodinger/pymol-open-source/blob/master/layer1/Color.cpp#L37
_color_cycle = [
    26,   # carbon
    5,    # cyan
    154,  # lightmagenta
    6,    # yellow
    9,    # salmon
    29,   # hydrogen
    11,   # slate
    13,   # orange
    10,   # lime
    5262, # deepteal
    12,   # hotpink
    36,   # yelloworange
    5271, # violetpurple
    124,  # grey70
    17,   # marine
    18,   # olive
    5270, # smudge
    20,   # teal
    5272, # dirtyviolet
    52,   # wheat
    5258, # deepsalmon
    5274, # lightpink
    5257, # aquamarine
    5256, # paleyellow
    15,   # limegreen
    5277, # skyblue
    5279, # warmpink
    5276, # limon
    53,   # violet
    5278, # bluewhite
    5275, # greencyan
    5269, # sand
    22,   # forest
    5266, # lightteal
    5280, # darksalmon
    5267, # splitpea
    5268, # raspberry
    104,  # grey50
    23,   # deepblue
    51,   # brown
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
        cmd.color(color, f"(e. C & c. {chain} & ({selection}))")
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
def cbattr(selection="(all)"):
    """
    DESCRIPTION
        Color by residue by attribute.
    USAGE
        cbattr [ selection ]
    """
    negative = ["ASP", "GLU"]
    positive = ["ARG", "HIS", "LYS"]
    hydrophobic = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
    hydrophilic = ["ASN", "GLN", "SER", "THR"]
    special = ["CYS", "GLY", "PRO"]
    cmd.color("red", f"({selection}) & resn {'+'.join(negative)}")
    cmd.color("marine", f"({selection}) & resn {'+'.join(positive)}")
    cmd.color("orange", f"({selection}) & resn {'+'.join(hydrophobic)}")
    cmd.color("green", f"({selection}) & resn {'+'.join(hydrophilic)}")
    cmd.color("gray", f"({selection}) & resn {'+'.join(special)}")
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
