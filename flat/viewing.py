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

from pymol import cmd, cgo
import numpy as np


@cmd.extend
def beautify(selection="all", mode=0, *, _self=cmd):
    """
    DESCRIPTION
        Display molecules in pretty representation.
    USAGE
        beautify [ selection [, mode ]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        mode = int: Display mode. {default: 0}
    """
    mode = int(mode)
    _self.hide("everything", selection)
    _self.show("wire", f"{selection} and ! polymer.protein")
    _self.show("cartoon", f"{selection} and polymer.protein")
    _self.show(
        ["licorice", "line"][mode],
        f"((byres ({selection})) & (sc. | (n. CA | n. N & r. PRO)))"
    )
    _self.hide(f"({selection} and hydro and (e. C extend 1))")
    _self.color("atomic", f"({selection}) and ! e. C")
    return


@cmd.extend
def filter(mode=1, *, _self=cmd):
    """
    DESCRIPTION
        Apply different (instagram) filter to molecules for rendering.
    USAGE
        filter [ mode ]
    ARGUMENTS
        mode = int: ... {default: 1}
    """
    mode = int(mode)
    if mode == 1:
        # Simple FLAT style
        # See: https://bionerdnotes.wordpress.com/2018/11/12/getting-high-quality-pictures-in-pymol/
        # cmd.color("cmky")
        cmd.set("alignment_as_cylinders", 1)
        cmd.set("cartoon_highlight_color", "gray75")
        cmd.set("cartoon_sampling", 14)
        cmd.set("depth_cue", 1)
        cmd.set("dot_as_spheres", 1)
        cmd.set("hash_max", 300)
        cmd.set("line_as_cylinders", 1)
        cmd.set("mesh_as_cylinders", 1)
        cmd.set("nb_spheres_quality", 3)
        cmd.set("nonbonded_as_cylinders", 1)
        cmd.set("ray_shadows", 0)
        cmd.set("ray_trace_color", "black")
        cmd.set("ray_trace_fog", 0.0)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ribbon_as_cylinders", 1)
        cmd.set("ribbon_sampling", 10)
        cmd.set("specular", 0.0)
        cmd.set("surface_quality", 1)
    elif mode == 2:
        # Pretty figure style
        cmd.bg_color("white")
        cmd.set("specular", 0)
        cmd.set("ambient", 0.3)
        cmd.set("direct", 1.0)
        cmd.set("ray_trace_mode", 1)
        cmd.set("stick_radius", 0.2)
        cmd.set("cartoon_tube_radius", 0.2)
        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_cylindrical_helices", 0)
        cmd.set("cartoon_flat_sheets", 1.0)
        cmd.set("cartoon_smooth_loops", 0)
        cmd.set("cartoon_highlight_color", "grey50")
        cmd.set("antialias", 1)
        cmd.set("hash_max", 300)
    else:
        # Set back to defaults
        cmd.reinitialize("settings")
    return


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
    """
    cmd.set("auto_zoom", 0)
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

    return


@cmd.extend
def bounding_box(selection="all", state=0, vis=1, color="yellow", quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Draws a bounding box around a given selection. 
    USAGE
        bounding_box [ selection [, state [, vis [, color [, quiet ]]]]]
    """
    state, vis, quiet = int(state), int(vis), int(quiet)
    rgb = _self.get_color_tuple(color)

    xyz = np.array(_self.get_coords(selection, state)).T

    # For what this vector magic is see reference bellow:
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
        print(" Bounding Box:", result, "A")

    if vis:
        linewidth = 2.00
        obj = [
            cgo.LINEWIDTH, linewidth,
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

    return  result


@cmd.extend
def gridbox(selection="(all)", grid=(5, 5, 5), color="white", lw=2.0, *, _self=cmd):
    """
    DESCRIPTION
        Given selection, draw a grid box around it.
    USAGE
        gridbox [selection [, grid [, color ]]] 
    ARGUMENTS
        selection = str: Atom selection {default: all}
        grid = list: Number of grids on X, Y, and Z axis {default: [5, 5, 5]}
        color = str: Color of the grid {default: white}
    """
    vmin, vmax = np.array(_self.get_extent(selection))

    if _self.is_string(grid):
        grid = _self.safe_list_eval(grid)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    grid = np.array(grid, dtype=int)
    if grid.size != 3:
        raise ValueError("Grid must be length 3")   

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

    name = _self.get_unused_name("gridbox")
    _self.load_cgo(obj, name, zoom=0)

    return name


# Autocomplete
cmd.auto_arg[0].update({
    "beautify": cmd.auto_arg[0]["zoom"],
    "bounding_box": cmd.auto_arg[0]["zoom"],
    "gridbox": cmd.auto_arg[0]["zoom"],
})
