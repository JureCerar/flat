# Copyright (C) 2023-2025 Jure Cerar
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

from pymol import cmd, cgo
from .shapes import *
import numpy as np

# TODO:
# - Converter from PDB names to GLYCAM (and vice versa)

_GLYCAN_COLORS = {
    "BOND": "gray",  # bond color
    "NAG": "flat.blue",
    "MAN": "flat.green",
    "BMA": "flat.green",
    "GAL": "flat.yellow",
    "FUC": "flat.red",
    "SIA": "flat.purple",
}
"""Color definitions for glycans"""

_GLYCAN_SHAPES = {
    "NAG": "cube",
    "MAN": "sphere",
    "BMA": "sphere",
    "GAL": "sphere",
    "FUC": "cone",
    "SIA": "octahedron",
}
"""Shape definitions for glycans"""


def _get_ring(selection):
    """Determine which atoms define the sugar rings"""
    ring = []
    # identify the oxygens that belong to the ring
    cmd.iterate(
        f"{selection} and n. C1 extend 1 and ({selection} and n. O* and ! n. O1*)",
        "ring.append(name)",
        space=locals()
    )
    if not ring:
        # NOTE: Can fail if ring does not start with C1 e.g. SIA
        cmd.iterate(
            f"{selection} and n. C2 extend 1 and ({selection} and n. O* and ! n. O1*)",
            "ring.append(name)",
            space=locals()
        )
    tmp = cmd.get_unused_name("tmp")
    try:
        for carbon in range(1, 10):
            num = cmd.select(
                tmp, f"not hydrogen and (neighbor ({selection} and n. C{carbon}))")
            if num > 2:
                ring.append(f"C{carbon}")
        while True:
            num = cmd.select(
                tmp, f"{selection} and n. {ring[0]} extend 1 and n. {ring[-1]}")
            if num == 0:
                ring.pop()
            else:
                break
    finally:
        cmd.delete(tmp)
    return ring


@cmd.extend
def snfg(selection="all", state=0, size=2.0, name=None, *, _self=cmd):
    """
    DESCRIPTION
        Show Standard Notation For Glycans (SNFG)
    USAGE
        snfg [ selection [, state [, size, [, name ]]]]
    ARGUMENTS
        selection = str: Atom selection {default: all}
        state = int: Object state (0 for current state) {default: 0}
        size = float: Size of output shapes {default: 2.0}
        name = str: Name of output group {default: "SNFG"}
    REFERENCE
        https://www.ncbi.nlm.nih.gov/glycans/snfg.html
    SEE ALSO
        _GLYCAN_COLORS, _GLYCAN_SHAPES
    """
    state, size = int(state), float(size)

    # Get residue list
    residues = []
    cmd.iterate(
        f"{selection} & name C1",
        "residues.append((model, segi, chain, resn, resi))",
        space=locals()
    )

    # Get bonds list
    bonds = list()
    for res in residues:
        neighbors = []
        cmd.iterate(
            "neighbor /{}/{}/{}/{}`{}".format(*res),
            "neighbors.append((model, segi, chain, resn, resi))",
            space=locals(),
        )
        for neighbor in neighbors:
            bond = sorted([res, neighbor])
            bonds.append(tuple(bond))

    # Calculate centers of each residue
    means, normals = dict(), dict()
    for res in residues:
        ring = _get_ring("/{}/{}/{}/{}`{}".format(*res))
        names = "+".join(ring)
        coord = cmd.get_coords(
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
        norm = np.mean(normal, axis=0)
        means[res] = mean
        normals[res] = norm

    obj = []

    # Draw bonds
    color = _GLYCAN_COLORS.get("BOND", "gray")
    color = cmd.get_color_tuple(color)
    for a1, a2 in set(bonds):
        # Ignore amino acid residues
        if any(_ not in residues for _ in (a1, a2)):
            continue
        try:
            xyz1, xyz2 = means[a1], means[a2]
            o = [cgo.CYLINDER, *xyz1, *xyz2, 0.5, *color, *color]
            obj.extend(o)
        except KeyError as e:
            print(f"Error: Unknown residue: {e}")

    # Draw shapes
    for res in residues:
        try:
            mean = means[res]
            norm = normals[res]
            resn = res[3]
            shape = _GLYCAN_SHAPES.get(resn, None)
            color = _GLYCAN_COLORS.get(resn, "white")
            # NOTE: 1.73 scale factor is there for all shapes to appear same size
            if shape == "cube":
                length = 1.73 * np.array([size, size, size])
                o = cube(mean, norm, 0, length, color)
            elif shape == "cone":
                o = cone(mean, norm, size, 1.73 * size, color)
            elif shape == "octahedron":
                o = octahedron(mean, norm, 0, 1.73 * size, color)
            else:
                o = sphere(mean, norm, size, color)
            obj.extend(o)
        except KeyError as e:
            print(f"Error: Unknown residue: {e}")

    if not name:
        name = cmd.get_unused_name("SNFG")
    cmd.load_cgo(obj, name, state=state, zoom=0)


# Autocomplete
cmd.auto_arg[0].update({
    "snfg": cmd.auto_arg[0]["zoom"],
})
