# Copyright (C) 2023-2024 Jure Cerar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from pymol import cmd, menu
import itertools

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


def hex2rgb(hex):
    """hex to rgb"""
    return tuple(int(hex[i:i+2], 16) for i in (2, 4, 6))


def rgb2hex(rgb):
    """rgb to hex"""
    return "0x%02x%02x%02x" % rgb


def mix(rgb1, rgb2, x=0):
    """Mix two colors"""
    return tuple(
        int((1 - x) * a + x * b) for a, b in zip(rgb1, rgb2)
    ) 


@cmd.extend
def cbo(selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Color by objects
    USAGE
        cbc [ selection ]
    """
    color = itertools.cycle(_color_cycle)
    for obj in cmd.get_names("public_objects", 0, selection):
        _self.color(next(color), f"({selection}) & %{obj} & e. C")
    return


@cmd.extend
def cbc(selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Color by chains
    USAGE
        cbc [ selection ]
    """
    colors = itertools.cycle(_color_cycle)
    for chain in _self.get_chains(selection):
        _self.color(next(colors), f"({selection}) & c. {chain} & e. C")
    return


@cmd.extend
def cbe(selection="(all)", *, _self=cmd, **kwargs):
    """
    DESCRIPTION
        Color by elements
    USAGE
        cbe [ selection [, <element>=<color> ]]
    EXAMPLE
        cbe (all), C=green, oxygen=red
    """
    _self.color("white", f"({selection}) & e. H")
    _self.color("blue", f"({selection}) & e. N")
    _self.color("red", f"({selection}) & e. O")
    _self.color("yellow", f"({selection}) & e. S")
    for element, color in kwargs.items():
        _self.color(color, f"(({selection}) & e. {element})")
    return


@cmd.extend
def cbattr(selection="(all)", *args, _self=cmd):
    """
    DESCRIPTION
        Color by residue attribute.
    USAGE
        cbattr [ selection [, attribute ]]
    NOTES
        Attributes: negative, positive, polar, nonpolar, and special
    """

    if "negative" in args or not args:
        res = ["ASP", "GLU"]
        _self.color("red", f"({selection}) & resn {'+'.join(res)}")

    if "positive" in args or not args:
        res = ["ARG", "HIS", "LYS"]
        _self.color("marine", f"({selection}) & resn {'+'.join(res)}")

    if "polar" in args or not args:
        res = ["ASN", "GLN", "SER", "THR"]
        _self.color("tv_green", f"({selection}) & resn {'+'.join(res)}")

    if "nonpolar" in args or not args:
        res = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
        _self.color("yellow", f"({selection}) & resn {'+'.join(res)}")

    if "special" in args or not args:
        res = ["CYS", "GLY", "PRO"]
        _self.color("gray", f"({selection}) & resn {'+'.join(res)}")

    return


@cmd.extend
def set_colors(scheme, *, _self=cmd):
    """
    DESCRIPTION
        Set elements color scheme.
    USAGE
        set_colors scheme
    SCHEMES
        jmol, [more coming soon]
    """
    if scheme == "jmol":
        _self.set("auto_color", 0)
        _self.set_color("hydrogen", [1.000, 1.000, 1.000])
        _self.set_color("carbon", [0.567, 0.567, 0.567])
        _self.set_color("nitrogen", [0.189, 0.315, 0.976])
        _self.set_color("oxygen", [1.000, 0.051, 0.051])
        _self.set_color("fluorine", [0.567, 0.882, 0.314])
        _self.set_color("sulfur", [1.000, 1.000, 0.189])
        _self.recolor()


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


# Color palettes for PyMOL API only
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
cmd.set_color("flat.purple", [162, 108, 252])
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

# Novartis NEW colors and palettes 
palette.nvs_blue = _Colors("nvs.novartis", "nvs.turquoise")
palette.nvs_blue_rev = _Colors("nvs.turquoise", "nvs.novartis")
palette.nvs_purple = _Colors("nvs.purple", "nvs.coral")
palette.nvs_purple_rev = _Colors("nvs.coral", "nvs.purple")
palette.nvs_coral = _Colors("nvs.coral", "nvs.amber")
palette.nvs_coral_rev = _Colors("nvs.amber", "nvs.coral")
cmd.set_color("nvs.novartis", [4, 96, 169])  #0460a9
cmd.set_color("nvs.carmine", [141, 21, 27])
cmd.set_color("nvs.sienna", [231, 74, 33])
cmd.set_color("nvs.apricot", [236, 154, 30])
cmd.set_color("nvs.deep", [0, 32, 104]) #002068
cmd.set_color("nvs.coral", [255, 88, 93])  # ff585d
cmd.set_color("nvs.purple", [143, 45, 222])  # 8f2dde
cmd.set_color("nvs.turquoise", [80, 226, 208])  # 50e2d0
cmd.set_color("nvs.amber", [255, 193, 0])  # ffc100
menu.all_colors_list.append(
    ("novartis", [
        ("036", "nvs.novartis"),
        ("501", "nvs.carmine"),
        ("921", "nvs.sienna"),
        ("961", "nvs.apricot"),
        ("014", "nvs.deep"),
        ("933", "nvs.coral"),
        ("518", "nvs.purple"),
        ("388", "nvs.turquoise"),
        ("970", "nvs.amber"),
    ]),
)

# AlphaFold colors and palettes
palette.alphafold = _Colors("af.high", "af.medium", "af.low", "af.crit")
palette.alphafold_rev = _Colors("af.crit", "af.low", "af.medium", "af.high")
cmd.set_color("af.high", [0, 83, 214])
cmd.set_color("af.medium", [101, 203, 243])
cmd.set_color("af.low", [255, 219, 19])
cmd.set_color("af.crit", [255, 125, 69])
menu.all_colors_list.append(
    ("alphafold", [
        ("038", "af.high"),
        ("379", "af.medium"),
        ("980", "af.low"),
        ("942", "af.crit"),
    ]),
)
