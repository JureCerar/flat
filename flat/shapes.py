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

from pymol import cgo, cmd, CmdException
from pymol.chempy import cpv


def cube(center=(0, 0, 0), normal=(0, 0, 1), length=(1, 1, 1),
         color="", *, _self=cmd):
    """
    DESCRIPTION
        Return a CGO rectangular cuboid object.
    USAGE 
        cube [ center [, normal [, length, [, color ]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        length = float3: object's dimensions {default: 1, 1, 1}
        color = string: object color {default: None}
    REFERENCE
        https://wiki.pymol.org/index.php/Cubes
    """
    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    def add_normal(xyz):
        return obj.extend(
            [cgo.NORMAL] + cpv.transform(matrix, xyz)
        )

    def add_vertex(xyz):
        return obj.extend(
            [cgo.VERTEX] + cpv.add(center, cpv.transform(matrix, xyz))
        )

    x, y, z = length
    for s in [-1, 1]:
        # X-axis
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLE_STRIP)
        if color:
            obj.extend([cgo.COLOR, *color])
        add_normal((s, 0, 0))
        add_vertex(cpv.scale((x,  y, -z), 0.5 * s))
        add_vertex(cpv.scale((x, -y, -z), 0.5 * s))
        add_vertex(cpv.scale((x,  y,  z), 0.5 * s))
        add_vertex(cpv.scale((x, -y,  z), 0.5 * s))
        obj.append(cgo.END)
        # Y-axis
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLE_STRIP)
        if color:
            obj.extend([cgo.COLOR, *color])
        add_normal((0, s, 0))
        add_vertex(cpv.scale((x, y, -z), 0.5 * s))
        add_vertex(cpv.scale((-x, y, -z), 0.5 * s))
        add_vertex(cpv.scale((x, y,  z), 0.5 * s))
        add_vertex(cpv.scale((-x, y,  z), 0.5 * s))
        obj.append(cgo.END)
        # Z-axis
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLE_STRIP)
        if color:
            obj.extend([cgo.COLOR, *color])
        add_normal((0, 0, s))
        add_vertex(cpv.scale((x,  y, z), 0.5 * s))
        add_vertex(cpv.scale((x, -y, z), 0.5 * s))
        add_vertex(cpv.scale((-x,  y, z), 0.5 * s))
        add_vertex(cpv.scale((-x, -y, z), 0.5 * s))
        obj.append(cgo.END)

    return obj


def sphere(center=(0, 0, 0), normal=(0, 0, 1), radius=1,
           color="", *, _self=cmd):
    """
    DESCRIPTION
        Return a CGO sphere object.
    USAGE 
        sphere [ center [, normal [, radius, [, color ]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        length = float: object's radius {default: 1}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    if color:
        obj.extend([cgo.COLOR, *color])
    obj.extend([cgo.SPHERE, *center, radius * 0.5])

    return obj


def cylinder(center=(0, 0, 0), normal=(0, 0, 1), radius=1,
             height=1, color="", *, _self=cmd):
    """
    DESCRIPTION
        Return a CGO cylinder object.
    USAGE 
        sphere [ center [, normal [, radius, [, height [, color ]]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        radius = float: object's radius {default: 1}
        height = float: object's height {default: 1}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    xyz1 = cpv.add(center, cpv.transform(matrix, [0, 0,  0.5*height]))
    xyz2 = cpv.add(center, cpv.transform(matrix, [0, 0, -0.5*height]))
    obj.extend([cgo.CYLINDER, *xyz1, *xyz2, radius * 0.5, *color, *color])

    return obj


def cone(center=(0, 0, 0), normal=(0, 0, 1), radius=1,
         height=1, color="", *, _self=cmd):
    """
    DESCRIPTION
        Return a CGO cone object.
    USAGE 
        sphere [ center [, normal [, radius, [, height [, color ]]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        radius = float: object's radius {default: 1}
        height = float: object's height {default: 1}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    xyz1 = cpv.add(center, cpv.transform(matrix, [0, 0,  0.5*height]))
    xyz2 = cpv.add(center, cpv.transform(matrix, [0, 0, -0.5*height]))
    obj.extend([cgo.CONE, *xyz1, *xyz2, radius * 0.5, 0, *color, *color, 1, 1])

    return obj


def polygon(center=(0, 0, 0), normal=(0, 0, 1), radius=1,
            sides=3, color="",  *, _self=cmd):
    """
    DESCRIPTION
        Return a CGO polygon object.
    """
    from math import sin, cos, pi

    if int(sides) < 3:
        raise CmdException()

    if color and isinstance(color, str):
        color = list(_self.get_color_tuple(color))
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    def add_normal(xyz):
        return obj.extend(
            [cgo.NORMAL] + cpv.transform(matrix, xyz)
        )

    def add_vertex(xyz):
        return obj.extend(
            [cgo.VERTEX] + cpv.add(center, cpv.transform(matrix, xyz))
        )

    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_FAN)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, 1))
    add_vertex((0, 0, 0))
    for i in range(sides + 1):
        angle = 2 * pi * i / sides
        x1, y1 = cos(angle) * radius, sin(angle) * radius
        add_vertex((x1, y1, 0))
    obj.append(cgo.END)

    return obj
