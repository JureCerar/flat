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
from typing import Iterable
import collections
import itertools
import colorsys
import re


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

# Attributes for coloring
_attributes = ["negative", "positive", "polar", "nonpolar", "special"]


# Regular expression for HEX colors
_HEX_COLOR_LONG = re.compile(r'^[0-9a-fA-F]{6}$')
_HEX_COLOR_SHORT = re.compile(r'^[0-9a-fA-F]{3}$')


class Color:
    """Converts and manipulates common PyMOL color representation: colorname, RGB, or HEX. PyMOL API only."""

    def __init__(self, color=None) -> None:
        self._rgb = self._transform(color)

    @staticmethod
    def _transform(color: None | str | tuple = None, *, _self=cmd) -> tuple[float, float, float]:
        """Transfrom input color name or color HEX to RGB float tuple."""
        if not color:
            rgb = (0., 0., 0.)

        elif isinstance(color, Color):
            rgb = color._rgb

        elif isinstance(color, str):
            color = color.lstrip('#')
            if _HEX_COLOR_LONG.match(color):
                rgb = tuple(int(color[i:i+2], 16) / 255 for i in (0, 2, 4))
            elif _HEX_COLOR_SHORT.match(color):
                rgb = tuple(int(color[i:i+1] * 2, 16) / 255 for i in (0, 1, 2))
            else:
                rgb = _self.get_color_tuple(color)
                if not rgb:
                    raise ValueError("Invalid color input")

        elif isinstance(color, Iterable):
            if len(color) != 3:
                raise ValueError("RGB color must contain 3 values")
            if all(isinstance(v, int) for v in color):
                if not all([0 <= v <= 255 for v in color]):
                    raise ValueError("RGB values must be between [0, 255]")
                rgb = tuple(v / 255 for v in color)
            elif all(isinstance(v, float) for v in color):
                if not all([0. <= v <= 1. for v in color]):
                    raise ValueError("RGB values must be between [0, 1]")
                rgb = tuple(color)
            else:
                raise ValueError("RGB values must be integer of float")

        else:
            return ValueError("Invalid color input")

        return rgb

    def __str__(self):
        return str(self._rgb)

    def __repr__(self):
        return repr(self._rgb)

    def __iter__(self):
        for v in self._rgb:
            yield v

    def __lt__(self, other):
        return self._rgb < other

    def __add__(self, other):
        if isinstance(other, (int, float)):
            color = tuple(min(max(0., v + other), 1.) for v in self)
        elif isinstance(other, Color):
            color = tuple(min(max(0., p + q), 1.) for p, q in zip(self, other))
        else:
            raise TypeError(
                "Unsupported operation for type: " + type(other).__name__)
        return Color(color)

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            color = tuple(min(max(0., v - other), 1.) for v in self)
        elif isinstance(other, Color):
            color = tuple(min(max(0., p - q), 1.) for p, q in zip(self, other))
        else:
            raise TypeError(
                "Unsupported operation for type: " + type(other).__name__)
        return Color(color)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            color = tuple(min(max(0, v * other), 1) for v in self)
        else:
            raise TypeError(
                "Unsupported operation for type: " + type(other).__name__)
        return Color(color)

    @staticmethod
    def rgb(r: int | float, g: int | float, b: int | float):
        """RGB color constructor method 
        Args:
            r (int, float): Red color fraction between [0, 255]
            g (int, float): Green color fraction between [0, 255]
            b (int, float): Blue color fraction between [0, 255]
        Returns:
            Color: Color class.
        Example:
            >>> Color.rgb(100, 200, 176)
            >>> Color.rgb(0.39, 0.78, 0.69)
        """
        return Color((r, g, b))

    @staticmethod
    def hex(string: str):
        """HEX color constructor method
        Args:
            string (str): Hex color string [0, 255]
        Returns:
            Color: Color class.
        Example:
            >>> Color.hex("64C8B0")
            >>> Color.hex("3F7")
        """
        return Color(string)

    def invert(self) -> None:
        """Inverts color
        Returns:
            None
        Examples:
            >>> Color.rgb(100, 200, 176).invert()
            (155, 55, 79)
        """
        self._rgb = tuple(1. - v for v in self)

    def to_rgb(self) -> tuple[int, int, int]:
        """Return color in RGB format
        Returns:
            tuple: List of HSL values.
        Example:
            >>> Color.hex('64c8b0').to_rgb()
            (100, 200, 176)
        """
        return tuple(round(v * 255) for v in self)

    def to_hex(self) -> str:
        """Return color in HEX format
        Returns:
            str: Color HEX string.
        Example:
            >>> Color.rgb(100, 200, 176).to_hex()
            '0x64c8b0'
        """
        return "0x%02x%02x%02x" % tuple(round(v * 255) for v in self)

    def to_hsl(self) -> tuple[float, float, float]:
        """Return color in HSL format
        Returns:
            tuple: List of HSL values.
        Example:
            >>> Color.rgb(100, 200, 176).to_hsl()
            (0.4600, 0.4762, 0.5882)
        """
        h, l, s = colorsys.rgb_to_hls(*self)
        return (h, s, l)


class Palette(list):
    """Class for better color palette handling. Behaves same as Python list. PyMOL API only."""

    def __init__(self, *args) -> None:
        values = (Color(v) for v in args)
        super().__init__(values)

    def sort(self, *, key: None = None, reverse: bool = False) -> None:
        """Sort the list in ascending order and return None.

        The sort is in-place (i.e. the list itself is modified) and stable (i.e. the
        order of two equal elements is maintained).

        If a key function is given, apply it once to each list item and sort them,
        ascending or descending, according to their function values.

        The reverse flag can be set to sort in descending order.

        Read more about [color sorting](https://www.alanzucconi.com/2015/09/30/colour-sorting).
        
        Args:
            key (expr): A function to specify the sorting criteria(s)
            reverse (bool): Will sort the list in descending order
        Example:
            >>> plt = Palette('0F0',  'F00', 'F0F', '00F', 'FF0', '0FF')
            >>> plt.sort(key=lambda x: colorsys.rgb_to_hls(*x))
            >>> print([c.to_hex() for c in plt])
            ['ff0000', 'ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff']    
        """
        super().sort(key=key, reverse=reverse)

    def append(self, object: str | tuple) -> None:
        """Append object to the end of the list"""
        super().append(Color(object))

    def extend(self, iterable: Iterable) -> None:
        """Extend list by appending elements from the iterable"""
        values = (Color(v) for v in iterable)
        super().extend(values)

    def insert(self, index: int, item: str | tuple) -> None:
        """Insert object before index"""
        super().insert(index, Color(item))

    def __setitem__(self, index, value):
        if isinstance(index, slice):
            values = [Color(v) for v in value]
            super().__setitem__(index, values)
        else:
            super().__setitem__(index, Color(value))

    def interp(self, frac: float) -> Color:
        """Computes the linear interpolation between all colors in palette, 
        if the parameter `frac` is inside [0, 1).
        Args:
            frac (float): Value between [0, 1] at which to evaluate the interpolated values color 
        Returns:
            result (Color): Interpolated Color object
        Examples:
            >>> plt = Palette('F00', '0F0', '00F')
            >>> plt.interp(0.25).to_hex()
            '808000'
            >>> nsteps = 10
            >>> for i in range(nsteps + 1): 
            >>>     plt.interp(i / nsteps).to_hex()
        """
        x = min(max(0., frac), 1.) * (len(self) - 1)
        if x.is_integer():
            value = self[int(x)]
        else:
            n, f = int(x), x - int(x)
            ca, cb = self[n], self[n+1]
            value = tuple((1 - f) * a + f * b for a, b in zip(ca, cb))
        return Color(value)



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
        _self.color(next(color), f"({selection}) & %{obj}")
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
        _self.color(next(colors), f"({selection}) & c. '{chain}'")
    return


@cmd.extend
def cbe(selection="(all)", *, _self=cmd, **kwargs):
    """
    DESCRIPTION
        Color by elements. By default skips coloring carbons.
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
        _self.color(color, f"({selection}) & e. {element}")
    return


@cmd.extend
def cbattr(selection="(all)", attribute=[], _self=cmd):
    """
    DESCRIPTION
        Color by residue attribute.
    USAGE
        cbattr [ selection [, attribute ]]
    ARGUMENTS
        selection = str: Atom selection {default: all}
        attribute = str|list: Attributes to color {default: None}
    NOTES
        Available attributes: 
            negative, positive, polar, nonpolar, and special
    """
    if not attribute:
        attribute = _attributes

    if _self.is_string(attribute):
        try:
            attribute = _self.safe_list_eval(attribute)
        except:
            attribute = attribute.split()

    if "negative" in attribute:
        res = ["ASP", "GLU"]
        _self.color("red", f"({selection}) & resn {'+'.join(res)}")
    if "positive" in attribute:
        res = ["ARG", "HIS", "LYS"]
        _self.color("marine", f"({selection}) & resn {'+'.join(res)}")
    if "polar" in attribute:
        res = ["ASN", "GLN", "SER", "THR"]
        _self.color("tv_green", f"({selection}) & resn {'+'.join(res)}")
    if "nonpolar" in attribute:
        res = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL"]
        _self.color("yellow", f"({selection}) & resn {'+'.join(res)}")
    if "special" in attribute:
        res = ["CYS", "GLY", "PRO"]
        _self.color("gray", f"({selection}) & resn {'+'.join(res)}")


@cmd.extend
def cbalpha(selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Color by C alpha atom.
    USAGE
        cbalpha [selection]
    """
    color_list = collections.defaultdict(lambda: 25) # gray is default
    _self.iterate(f"bca. ({selection})",
                  "color_list[object,segi,chain,resi] = color", space=locals())
    _self.alter(selection,
                "color = color_list[object,segi,chain,resi]", space=locals())
    _self.recolor()
    return


@cmd.extend
def cbs(selection="all", palette="rainbow", repr="ribbon", *, _self=cmd):
    """
    DESCRIPTION
        Color object by state
    USAGE
        cbs [ selection [, palette [, repr ]]]
    ARGUMENTS
        selection = str: Atom selection {default: all}
        palette = str|list: List of colors {default: rainbow}
        repr = str: Representation to color {default: ribbon}
    EXAMPLE
        >>> set all_states, 1
        >>> show_as ribbon
        >>> cbs all
        >>> cbs all, red blue
        >>> cbs all, repr=cartoon
        >>> cbs all, undo
    NOTES
        Due to how PyMol handles states only one representation can be
        colored at the time. To undo the coloring use `undo` as palette. 
    """
    num_states = _self.count_states(selection)

    # Special palette keywords:
    if _self.is_string(palette):
        if palette == "auto":
            # Do auto state coloring
            for i in range(num_states):
                _self.set(f"{repr}_color", "auto", selection, state=i+1)
            return
        elif palette == "undo":
            # Reset representation do default
            for i in range(num_states):
                _self.set(f"{repr}_color", -1, selection, state=i+1)
            return
        elif palette == "rainbow":
            palette = ["red", "orange", "yellow", "green", "blue", "purple"]
        elif palette == "rainbow_rev":
            palette = ["purple", "blue", "green", "yellow", "orange", "red"]
        else:
            try:
                palette = _self.safe_list_eval(palette)
            except:
                palette = palette.split()

    if len(palette) < 2:
        raise IndexError("Specify at least two colors")

    # Mix colors one per state
    color_hex = [_self.get_color_tuple(c) for c in palette]
    for i in range(num_states):
        x = i / (num_states - 1) * (len(palette) - 1)
        if x.is_integer():
            color = color_hex[int(x)]
        else:
            n, f = int(x), x - int(x)
            c1, c2 = color_hex[n], color_hex[n+1]
            color = tuple((1-f) * a + f * b for a, b in zip(c1, c2))
        _self.set(f"{repr}_color", color, selection, state=i+1)


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


@cmd.extend
def get_colors(selection="(all)", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Get color of atoms in selection. Returns list of indices and RGB tuples.
    USAGE
        get_colors [selection [, quiet ]]
    PYMOL API
        get_colors(selection) -> [(index, (r, g, b)), ... ]
    """
    color_list = []
    get_color_tuple = _self.get_color_tuple
    _self.iterate(selection,
                  "color_list.append((index, get_color_tuple(color)))", space=locals())
    if not quiet:
        for index, color in color_list:
            print(index, color)
    return color_list
    

# Autocompletion
cmd.auto_arg[0].update({
    "cbo": cmd.auto_arg[0]["zoom"],
    "cbc": cmd.auto_arg[0]["zoom"],
    "cbe": cmd.auto_arg[0]["zoom"],
    "cbattr": cmd.auto_arg[0]["zoom"],
    "cbalpha": cmd.auto_arg[0]["zoom"],
    "cbs": cmd.auto_arg[0]["zoom"],
    "get_colors": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "cbattr": [cmd.Shortcut(_attributes), "arg", ""],
    "cbs": cmd.auto_arg[0]["color"],
})
cmd.auto_arg[2].update({
    "cbs": cmd.auto_arg[0]["as"]
})


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

# Novartis colors and palettes 
cmd.set_color("nvs.blue", [4, 96, 169])  #0460a9
cmd.set_color("nvs.deep", [0, 32, 104])  #002068
cmd.set_color("nvs.coral", [255, 88, 93])  #ff585d
cmd.set_color("nvs.amber", [255, 193, 0])  #ffc100
cmd.set_color("nvs.light", [208, 208, 208])  #d0d0d0
cmd.set_color("nvs.gray", [167, 167, 167])  #a7a7a7
menu.all_colors_list.append(
    ("novartis", [
        ("036", "nvs.blue"),
        ("014", "nvs.deep"),
        ("933", "nvs.coral"),
        ("970", "nvs.amber"),
        ("777", "nvs.light"),
        ("555", "nvs.gray"),
    ]),
)

# AlphaFold colors and palettes
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
