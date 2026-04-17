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
:mod:`flat.viewing`
===================
Module to help with molecule visualization
"""

from pymol import cmd, CmdException, cgo
from pymol.chempy import cpv
import numpy as np


@cmd.extend
def beautify(selection="all", mode=0, *, _self=cmd):
    """
    DESCRIPTION
        Display molecules in pretty representation.
    USAGE
        beautify [ selection [, mode ]]
    ARGUMENTS
        selection : str, optional
            Atom selection.  
        mode : int, default = 0
            Display mode.
    """
    mode = int(mode)
    _self.hide("everything", selection)
    _self.show("wire", f"{selection} and ! polymer.protein")
    _self.show("cartoon", f"{selection} and polymer.protein")
    if mode == 1:
        repr = "licorice"
    elif mode == 2:
        repr = "wire"
    else:
        repr = "line"
    _self.show(repr, f"((byres ({selection})) & (sc. | (n. CA | n. N & r. PRO)))")
    _self.hide(f"({selection} and hydro and (e. C extend 1))")
    _self.color("atomic", f"({selection}) and ! e. C")


@cmd.extend
def title(label, selection="polymer", size=26, va=0.8, ha=0.0, *, _self=cmd):
    """
    DESCRIPTION
        Display title over selection.
    USAGE
        title label [ selection [, size [, va  [, ha ]]]]
    ARGUMENTS
        label : str
            Text to display
        selection : str, default = 'polymer'
            Atom selection
        size : int, default = 26
            Text size
        va : float, default = 0.0
            Vertical alignment where 0 is center
        ha : float, default = 0.8
            Horizontal alignment where 0 is center
    """
    size, ha, va = int(size), float(ha), float(va)
    com = _self.centerofmass(selection)
    vmin, vmax = _self.get_extent(selection)
    dist = cpv.distance(vmin, vmax) / 2
    name = _self.get_unused_name("title")
    _self.pseudoatom(name, label=label, pos=com)
    _self.set("label_position", (dist * ha, dist * va, 0), name)
    _self.set("label_size", size, name)


@cmd.extend
def filter(mode=1, *, _self=cmd):
    """
    DESCRIPTION
        Apply different preset filters to molecules for rendering:

        * 0 = PyMol default style 
        * 1 = Simple FLAT style  
        * 2 = Pretty figure style
    USAGE
        filter [ mode ]
    ARGUMENTS
        mode : int, default = 1
            Render mode preset
    """
    mode = int(mode)
    if mode == 1:
        # Simple FLAT style
        # See: https://bionerdnotes.wordpress.com/2018/11/12/getting-high-quality-pictures-in-pymol/
        # cmd.color("cmky")
        _self.set("alignment_as_cylinders", 1)
        _self.set("cartoon_highlight_color", "gray75")
        _self.set("cartoon_sampling", 14)
        _self.set("depth_cue", 1)
        _self.set("dot_as_spheres", 1)
        _self.set("hash_max", 300)
        _self.set("line_as_cylinders", 1)
        _self.set("mesh_as_cylinders", 1)
        _self.set("nb_spheres_quality", 3)
        _self.set("nonbonded_as_cylinders", 1)
        _self.set("ray_shadows", 0)
        _self.set("ray_trace_color", "black")
        _self.set("ray_trace_fog", 0.0)
        _self.set("ray_trace_mode", 1)
        _self.set("ribbon_as_cylinders", 1)
        _self.set("ribbon_sampling", 10)
        _self.set("specular", 0.0)
        _self.set("surface_quality", 1)
    elif mode == 2:
        # Pretty figure style
        _self.bg_color("white")
        _self.set("specular", 0)
        _self.set("ambient", 0.3)
        _self.set("direct", 1.0)
        _self.set("ray_trace_mode", 1)
        _self.set("stick_radius", 0.2)
        _self.set("cartoon_tube_radius", 0.2)
        _self.set("cartoon_fancy_helices", 1)
        _self.set("cartoon_cylindrical_helices", 0)
        _self.set("cartoon_flat_sheets", 1.0)
        _self.set("cartoon_smooth_loops", 0)
        _self.set("cartoon_highlight_color", "grey50")
        _self.set("antialias", 1)
        _self.set("hash_max", 300)
    else:
        # Set back to defaults
        _self.reinitialize("settings")


class PutCenterCallback(object):
    prev_v = None

    def __init__(self, name, corner=0):
        self.name = name
        self.corner = corner
        self.cb_name = cmd.get_unused_name("_cb")
        return

    def load(self):
        cmd.load_callback(self, self.cb_name)
        return

    def __call__(self):
        if self.name not in cmd.get_names("objects"):
            import threading
            threading.Thread(None, cmd.delete, args=(self.cb_name,)).start()
            return

        v = cmd.get_view()
        if v == self.prev_v:
            return
        self.prev_v = v
        t = v[12:15]

        if self.corner:
            vp = cmd.get_viewport()
            R_mc = [v[0:3], v[3:6], v[6:9]]
            off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
            if self.corner in [2, 3]:
                off_c[0] *= -1
            if self.corner in [3, 4]:
                off_c[1] *= -1

            off_m = [
                R_mc[0][0]*off_c[0] + R_mc[0][1] *
                off_c[1] + R_mc[0][2]*off_c[2],
                R_mc[1][0]*off_c[0] + R_mc[1][1] *
                off_c[1] + R_mc[1][2]*off_c[2],
                R_mc[2][0]*off_c[0] + R_mc[2][1] *
                off_c[1] + R_mc[2][2]*off_c[2],
            ]
            t = [
                t[0] + off_m[0],
                t[1] + off_m[1],
                t[2] + off_m[2],
            ]

        z = -v[11] / 30.0
        m = [
            z, 0, 0, 0, 0, z, 0, 0, 0, 0, z,
            0, t[0] / z, t[1] / z, t[2] / z, 1
        ]
        cmd.set_object_ttt(self.name, m)
        return


@cmd.extend
def axis(name="axis", *, _self=cmd):
    """
    DESCRIPTION
        Puts coordinate axes to the lower left corner of the viewport.
    USAGE
        axis [ name ]
    ARGUMENTS
        name : str, default = 'axis'
            Name of the axis object 
    """
    _self.set("auto_zoom", 0)
    width = 0.06  # cylinder width
    length = 0.75  # cylinder length
    hight = 0.25  # cone hight
    d = width * 1.618  # cone base diameter
    obj = [
        cgo.CYLINDER, 0.0, 0.0, 0.0, length, 0.0, 0.0, width, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0, length, 0.0, width, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
        cgo.CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0, length, width, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
        cgo.CONE, length, 0.0, 0.0, hight +
        length, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
        cgo.CONE, 0.0, length, 0.0, 0.0, hight +
        length, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        cgo.CONE, 0.0, 0.0,   length, 0.0, 0.0, hight +
        length, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0
    ]
    # Display
    PutCenterCallback(name, 1).load()
    _self.load_cgo(obj, name, zoom=0)


@cmd.extend
def bounding_box(selection="all", state=0, vis=1, color="yellow", lw=2.0, *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Draws a bounding box around a given selection. 
    USAGE
        bounding_box [ selection [, state [, vis [, color ]]]]
    ARGUMENTS
        selection : str, optional
            Atom selection
        state : int, default = 0
            Object state where 0 is all states
        vis : bool, default = True
            Whether to visualize the bounding box
        color : str, default = 'yellow'
            Color of the bounding box
    RETURNS
        : ndarray
            Array of shape (3,) with bounding box dimensions in Angstroms   
    """
    state, vis, lw, quiet = int(state), int(vis), float(lw), int(quiet)
    rgb = _self.get_color_tuple(color)

    xyz = np.array(_self.get_coords(selection, state)).T

    # For what all this vector magic is see the reference bellow:
    #  _notebooks/2021-04-20-3D-Oriented-Bounding-Box.ipynb
    means = np.mean(xyz, axis=1)
    cov = np.cov(xyz)
    eval, evec = np.linalg.eig(cov)

    centered_xyz = xyz - means[:, np.newaxis]
    aligned_xyz = np.matmul(evec.T, centered_xyz)

    xmin, xmax = np.min(aligned_xyz[0, :]), np.max(aligned_xyz[0, :])
    ymin, ymax = np.min(aligned_xyz[1, :]), np.max(aligned_xyz[1, :])
    zmin, zmax = np.min(aligned_xyz[2, :]), np.max(aligned_xyz[2, :])

    box_xyz = np.array([
        [xmin, xmin, xmax, xmax, xmin, xmin, xmax, xmax],
        [ymin, ymax, ymax, ymin, ymin, ymax, ymax, ymin],
        [zmin, zmin, zmin, zmin, zmax, zmax, zmax, zmax]
    ])

    realigned_coords = np.matmul(evec, aligned_xyz)
    realigned_coords += means[:, np.newaxis]

    box = np.matmul(evec, box_xyz)
    box += means[:, np.newaxis]

    result = [
        np.linalg.norm(box[:, 4] - box[:, 7]),
        np.linalg.norm(box[:, 3] - box[:, 7]),
        np.linalg.norm(box[:, 6] - box[:, 7]),
    ]

    if not quiet:
        print(" Bounding Box:", result, "Å")

    if vis:
        obj = [
            cgo.LINEWIDTH, lw,
            cgo.BEGIN, cgo.LINES,
            cgo.COLOR, *rgb,
            # z1 plane boundary
            cgo.VERTEX, *box[:, 0], cgo.VERTEX, *box[:, 1],
            cgo.VERTEX, *box[:, 1], cgo.VERTEX, *box[:, 2],
            cgo.VERTEX, *box[:, 2], cgo.VERTEX, *box[:, 3],
            cgo.VERTEX, *box[:, 3], cgo.VERTEX, *box[:, 0],
            # z2 plane boundary
            cgo.VERTEX, *box[:, 4], cgo.VERTEX, *box[:, 5],
            cgo.VERTEX, *box[:, 5], cgo.VERTEX, *box[:, 6],
            cgo.VERTEX, *box[:, 6], cgo.VERTEX, *box[:, 7],
            cgo.VERTEX, *box[:, 7], cgo.VERTEX, *box[:, 4],
            # z1 and z2 connecting boundaries
            cgo.VERTEX, *box[:, 0], cgo.VERTEX, *box[:, 4],
            cgo.VERTEX, *box[:, 1], cgo.VERTEX, *box[:, 5],
            cgo.VERTEX, *box[:, 2], cgo.VERTEX, *box[:, 6],
            cgo.VERTEX, *box[:, 3], cgo.VERTEX, *box[:, 7],
            cgo.END
        ]
        name = _self.get_unused_name("BoundingBox")
        _self.load_cgo(obj, name, zoom=0)

    return result


@cmd.extend
def gridbox(selection="all", grid=(5, 5, 5), color="yellow", lw=2.0, *, _self=cmd):
    """
    DESCRIPTION
        Draw a grid box around selection.
    USAGE
        gridbox [ selection [, grid [, color ]]]
    ARGUMENTS
        selection : str, optional
            Atom selection
        grid : array-like, shape (3,), default = (5, 5, 5)
            Number of grids on X, Y, and Z axis
        color : str, default = 'yellow'
            Color of the grid box
    """
    vmin, vmax = np.array(_self.get_extent(selection))

    if _self.is_string(grid):
        grid = _self.safe_list_eval(grid)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    grid = np.array(grid, dtype=int)
    if grid.size != 3:
        raise CmdException("Grid must be length 3")

    dv = (vmax - vmin) / grid
    obj = [
        cgo.LINEWIDTH, float(lw),
        cgo.BEGIN, cgo.LINES,
    ]
    if color:
        obj.extend([cgo.COLOR, *color])

    for i in range(grid[0]):
        for j in range(grid[1]):
            for k in range(grid[2]):
                obj.extend([
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+k*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+i*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+j*dv[1], vmin[2]+(k+1)*dv[2],
                    cgo.VERTEX, vmin[0]+(i+1)*dv[0], vmin[1]+(j+1)*dv[1], vmin[2]+(k+1)*dv[2],
                ])
    obj.append(cgo.END)
    # Draw object
    name = _self.get_unused_name("gridbox")
    _self.load_cgo(obj, name, zoom=0)


# Autocomplete
cmd.auto_arg[0].update({
    "beautify": cmd.auto_arg[0]["zoom"],
    "bounding_box": cmd.auto_arg[0]["zoom"],
    "gridbox": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "title": cmd.auto_arg[0]["zoom"],
})
