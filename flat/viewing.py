#!/usr/bin/env python3
from pymol import cmd, cgo
from chempy import cpv


@cmd.extend
def beautify(selection="all", mode=0):
    """
    DESCRIPTION
        Display molecules in pretty representation.
    USAGE
        beautify [ selection [, mode ]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        mode = int: ... {default: 0}
    """
    mode = int(mode)
    cmd.hide("everything")
    cmd.show("cartoon", f"{selection} and polymer.protein")
    cmd.show("lines", f"((byres ({selection})) & (sc. | (n. CA | n. N & r. PRO)))")
    cmd.hide(f"({selection} and hydro and (e. C extend 1))")
    cmd.show("wire", f"{selection} and ! polymer.protein")
    # cmd.util.cbc()
    cmd.color("atomic", f"({selection}) and not elem C")
    return None


@cmd.extend
def filter(mode=1):
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
            off_m = cpv.transform(R_mc, off_c)
            t = cpv.add(t, off_m)

        z = -v[11] / 30.0
        m = [
            z, 0, 0, 0, 0, z, 0, 0, 0, 0, z,
            0, t[0] / z, t[1] / z, t[2] / z, 1
        ]
        cmd.set_object_ttt(self.name, m)
        return


@cmd.extend
def axis(name="axis"):
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
    cmd.load_cgo(obj, name)
    return


# https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/Draw_Protein_Dimensions.py
def dimensions():
    return
