# Copyright (C) 2023-2026 Jure Cerar
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

"""
:mod:`flat.glycans`
===================
Module for visualizing and editing glycans.
"""

from pymol import cmd, cgo
from .shapes import *
import numpy as np


# See: https://doi.org/10.1107/S2059798316016910
_GLYCAN_ONE_LETTER = {
    "ARA": "A", "AHR": "A", "LDY": "D", "RIB": "R", "RIP": "R", "XYS": "X",
    "XYP": "X", "GLA": "L", "GAL": "L", "GLC": "G", "BGC": "G", "MAN": "M",
    "BMA": "M", "FRU": "C", "FUC": "F", "FUL": "F", "GTR": "O", "ADA": "O",
    "GCU": "Z", "BDP": "Z", "IDR": "U", "A2G": "V", "NGA": "V", "NDG": "Y",
    "NAG": "Y", "BM3": "W", "SIA": "S", "SLB": "S",
}


_GLYCAN_REPR = {
    "Y": ("cube", "flat.blue"),  # Man
    "M": ("sphere", "flat.green"),  # GlcNAc
    "V": ("cube", "flat.yellow"),  # GalNAc
    "L": ("sphere", "flat.yellow"),  # Gal
    "F": ("cone", "flat.red"),  # Fuc
    "S": ("octahedron", "flat.purple"),  # Sia
}


def _get_glycan(key: str):
    """Try to get one-letter code for glycan key"""
    key = key.strip()

    if len(key) == 1:
        return key
    
    elif len(key) == 3:
        if key in _GLYCAN_ONE_LETTER:
            return _GLYCAN_ONE_LETTER[key]
        # Check if GLYCAM notation
        if key[0] in list("012346ZYXWVUTSRQP"):
            if key[2] in list("AB"):
                return key[1]
                
    return None


def _get_ring(selection: str, _self=cmd):
    """Determine which atoms define the glycan rings"""
    ring = []
    # identify the oxygens that belong to the ring
    _self.iterate(
        f"{selection} and n. C1 extend 1 and ({selection} and n. O* and ! n. O1*)",
        "ring.append(name)",
        space=locals()
    )
    if not ring:
        # NOTE: Can fail if ring does not start with C1 e.g. SIA
        _self.iterate(
            f"{selection} and n. C2 extend 1 and ({selection} and n. O* and ! n. O1*)",
            "ring.append(name)",
            space=locals()
        )
    tmp = _self.get_unused_name("tmp")
    try:
        for carbon in range(1, 10):
            num = _self.select(
                tmp, f"not hydrogen and (neighbor ({selection} and n. C{carbon}))")
            if num > 2:
                ring.append(f"C{carbon}")
        while True:
            num = _self.select(
                tmp, f"{selection} and n. {ring[0]} extend 1 and n. {ring[-1]}")
            if num == 0:
                ring.pop()
            else:
                break
    finally:
        _self.delete(tmp)
    return ring


@cmd.extend
def snfg(selection="all", state=0, size=2.0, width=0.5, name="SNFG",  *, shape_color=None, _self=cmd):
    """
    DESCRIPTION
        Show Standard Notation For Glycans (SNFG).
    USAGE
        snfg [ selection [, state [, size, [, width [, name ]]]]]
    ARGUMENTS
        selection : str, optional
            Atom selection.
        state : int, default = 0
            Object state. 0 for current state.
        size : float, default = 2.0
            Size of output shapes.
        width : float, default = 0.5
            Width of bonds between shapes {default: 0.5}
        name : str, default = 'SNFG'
            Name of output group.
    REFERENCE
        https://www.ncbi.nlm.nih.gov/glycans/snfg.html
    """
    state, size, width = int(state), float(size), float(width)

    # Scale factor is there for all shapes to _appear_ the same size
    SCALE = 1.73

    # Extend glycan dictionary
    _glycan_repr = _GLYCAN_REPR.copy()
    if shape_color:
        _glycan_repr.update(shape_color)

    # Get residue list
    residues = []
    _self.iterate(
        f"{selection} & name C1",
        "residues.append((model, segi, chain, resn, resi))",
        space=locals()
    )

    # Get bonds list
    bonds = list()
    for res in residues:
        neighbors = []
        _self.iterate(
            "neighbor /{}/{}/{}/{}`{}".format(*res),
            "neighbors.append((model, segi, chain, resn, resi))",
            space=locals(),
        )
        for neighbor in neighbors:
            bond = sorted([res, neighbor])
            bonds.append(tuple(bond))

    # Process each glycan: center, normal, color, and shape
    means, normals = dict(), dict()
    colors, shapes = dict(), dict()
    for res in residues:
        ring = _get_ring("/{}/{}/{}/{}`{}".format(*res), _self)
        names = "+".join(ring)
        coord = _self.get_coords(
            "/{}/{}/{}/{}`{} & name {}".format(*res, names),
            state=state,
        )
        mean = np.mean(coord, axis=0)
        normal = []
        for i in range(len(coord) - 1):
            normal.append(0.5 * (
                np.cross((coord[i] - coord[i - 1]), (mean - coord[i])) +
                np.cross((coord[i + 1] - coord[i]), (mean - coord[i]))
            ))
        # NOTE: Invert normal
        means[res] = mean
        normals[res] = -np.mean(normal, axis=0)
        try:
            key = _get_glycan(res[3])
            shape, color = _glycan_repr.get(key)
        except KeyError as e:
            print(f"Error: Unknown residue: {e}")
            # Default settings for unknown residues
            shape, color = (None, "gray")
        colors[res] = _self.get_color_tuple(color)
        shapes[res] = shape
        
    # Objects for drawing
    obj = []

    # Draw bonds
    for a1, a2 in set(bonds):
        # Ignore amino acid residues
        if any(_ not in residues for _ in (a1, a2)):
            continue
        try:
            xyz1, xyz2 = means[a1], means[a2]
            color1, color2 = colors[a1], colors[a2]
            o = [cgo.CYLINDER, *xyz1, *xyz2, width, *color1, *color2]
            obj.extend(o)
        except KeyError as e:
            print(f"Error: Unknown residue: {e}")

    # Draw shapes
    for res in residues:
        try:
            mean = means[res]
            norm = normals[res]
            shape = shapes[res]
            color = colors[res]
            if shape == "cube":
                length = SCALE * np.array([size, size, size])
                o = cube(mean, norm, 0, length, color)
            elif shape == "cone":
                o = cone(mean, norm, size, SCALE * size, color)
            elif shape == "octahedron":
                o = octahedron(mean, norm, 0, SCALE * size, color)
            elif shape == "star":
                o = star(mean, norm, 0, 5, (size, 0.5 * size), SCALE * size, color)
            else:
                o = sphere(mean, norm, size, color)
            obj.extend(o)
        except KeyError as e:
            print(f"Error: Unknown residue: {e}")

    if not name:
        name = _self.get_unused_name("SNFG")
    _self.load_cgo(obj, name, state=state, zoom=0)


# Autocomplete
cmd.auto_arg[0].update({
    "snfg": cmd.auto_arg[0]["zoom"],
})
