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

from pymol import cgo, cmd
from pymol.chempy import cpv


@cmd.extend
def cube(center=(0, 0, 0), normal=(0, 0, 1), rotation=0, length=(1, 1, 1),
         color="", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO rectangular cuboid object.
    USAGE 
        cube [ center [, normal [, rotation, [, length [, color ]]]]]
    ARGUMENTS
        center = float3: shape position {default: 0, 0, 0}
        normal = float3: orientation of shape {default: 1, 0, 0}
        rotation = float: rotation around normal in rad {default: 0}
        length = float3: dimensions {default: 1, 1, 1}
        color = string: color name {default: None}
    REFERENCE
        https://wiki.pymol.org/index.php/Cubes
    """
    quiet = int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if _self.is_string(length):
        length = _self.safe_list_eval(length)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))
    rotmat = cpv.rotation_matrix(rotation, cpv.normalize(normal))

    def add_normal(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.NORMAL] + _xyz)

    def add_vertex(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.VERTEX] + cpv.add(center, _xyz))

    x, y, z = cpv.scale(length, 0.5)

    # Right
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((1, 0, 0))
    add_vertex(( x, -y,  z))
    add_vertex(( x, -y, -z))
    add_vertex(( x,  y,  z))
    add_vertex(( x,  y, -z))
    obj.append(cgo.END)

    # Left
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((-1, 0, 0))
    add_vertex((-x, -y,  z))
    add_vertex((-x,  y,  z))
    add_vertex((-x, -y, -z))
    add_vertex((-x,  y, -z))
    obj.append(cgo.END)

    # Top
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 1, 0))
    add_vertex((-x,  y,  z))
    add_vertex(( x,  y,  z))
    add_vertex((-x,  y, -z))
    add_vertex(( x,  y, -z))
    obj.append(cgo.END)

    # Bottom
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, -1, 0))
    add_vertex((-x, -y,  z))
    add_vertex((-x, -y, -z))
    add_vertex(( x, -y,  z))
    add_vertex(( x, -y, -z))
    obj.append(cgo.END)

    # Front
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, 1))
    add_vertex((-x, -y,  z))
    add_vertex(( x, -y,  z))
    add_vertex((-x,  y,  z))
    add_vertex(( x,  y,  z))
    obj.append(cgo.END)

    # Back
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_STRIP)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, -1))
    add_vertex((-x, -y, -z))
    add_vertex((-x,  y, -z))
    add_vertex(( x, -y, -z))
    add_vertex(( x,  y, -z))
    obj.append(cgo.END)

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def sphere(center=(0, 0, 0), normal=(0, 0, 1), radius=0.5,
           color="", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO sphere object.
    USAGE 
        sphere [ center [, normal [, radius, [, color ]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        radius = float: object's radius {default: 1}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    quiet = int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)

    obj = []

    if color:
        obj.extend([cgo.COLOR, *color])
    obj.extend([cgo.SPHERE, *center, radius])

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def cylinder(center=(0, 0, 0), normal=(0, 0, 1), radius=0.5,
             height=1.0, color="", *, quiet=1, _self=cmd):
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
    quiet = int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    else:
        color = (1, 1, 1)  #ffffff

    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    xyz1 = cpv.add(center, cpv.transform(matrix, [0, 0,  0.5 * height]))
    xyz2 = cpv.add(center, cpv.transform(matrix, [0, 0, -0.5 * height]))
    obj.extend([cgo.CYLINDER, *xyz1, *xyz2, radius, *color, *color])

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def cone(center=(0, 0, 0), normal=(0, 0, 1), radius=0.5,
         height=1.0, color="", *, quiet=1, _self=cmd):
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
    quiet = int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    else:
        color = (1, 1, 1)  #ffffff
    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    xyz1 = cpv.add(center, cpv.transform(matrix, [0, 0,  0.5*height]))
    xyz2 = cpv.add(center, cpv.transform(matrix, [0, 0, -0.5*height]))
    obj.extend([cgo.CONE, *xyz1, *xyz2, 0, radius, *color, *color, 1, 1])

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def polygon(center=(0, 0, 0), normal=(0, 0, 1), rotation=0, radius=0.5,
            sides=3, color="", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO n-sided polygon object.
    """
    from math import sin, cos, pi

    sides, quiet = int(sides), int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)

    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))
    rotmat = cpv.rotation_matrix(rotation, cpv.normalize(normal))

    def add_normal(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.NORMAL] + _xyz)

    def add_vertex(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.VERTEX] + cpv.add(center, _xyz))

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

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def prism(center=(0, 0, 0), normal=(0, 0, 1), rotation=0, radius=0.5,
          height=1, sides=3, color="",  *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO n-sided prism object.
    USAGE 
        sphere [ center [, normal [, rotation [, radius [, height [, sides  [, color ]]]]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        rotation = float: rotation around normal in rad {default: 0}
        radius = float: object's radius {default: 1}
        height = float: object's height {default: 1}
        sides = int: num of polygon sides {default: 3}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    from math import sin, cos, pi

    sides, quiet = int(sides), int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)

    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))
    rotmat = cpv.rotation_matrix(rotation, cpv.normalize(normal))

    def add_normal(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.NORMAL] + _xyz)

    def add_vertex(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.VERTEX] + cpv.add(center, _xyz))

    # Top
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_FAN)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, 1))
    add_vertex((0, 0, height * 0.5))
    for i in range(sides + 1):
        angle = 2 * pi * i / sides
        x1, y1 = cos(angle) * radius, sin(angle) * radius
        add_vertex((x1, y1, height * 0.5))
    obj.append(cgo.END)

    # Bottom
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_FAN)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, -1))
    add_vertex((0, 0, -height * 0.5))
    for i in range(sides + 1)[::-1]:
        angle = 2 * pi * i / sides
        x1, y1 = cos(angle) * radius, sin(angle) * radius
        add_vertex((x1, y1, -height * 0.5))
    obj.append(cgo.END)

    # Sides
    for i in range(sides + 1):
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLE_STRIP)
        if color:
            obj.extend([cgo.COLOR, *color])
        angle1 = 2 * pi * i / sides
        x1, y1 = cos(angle1) * radius, sin(angle1) * radius
        angle2 = 2 * pi * (i + 1) / sides
        x2, y2 = cos(angle2) * radius, sin(angle2) * radius
        add_normal(((x1 + x2) * 0.5, (y1 + y2) * 0.5, 0))
        add_vertex((x1, y1, -height * 0.5))
        add_vertex((x2, y2, -height * 0.5))
        add_vertex((x1, y1,  height * 0.5))
        add_vertex((x2, y2,  height * 0.5))
        obj.append(cgo.END)

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def pyramid(center=(0, 0, 0), normal=(0, 0, 1), rotation=0, radius=0.5,
            height=1.0, sides=3, color="",  *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO n-sided pyramid object.
    USAGE 
        pyramid [ center [, normal [, rotation, [, radius [, height [, sides  [, color ]]]]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: orientation of object  {default: 0, 0, 1}
        rotation = float: rotation around normal in rad {default: 0}
        radius = float: object's radius {default: 1}
        height = float: object's height {default: 1}
        sides = int: num of pyramid sides {default: 3}
        color = string: object color {default: None}
    REFERENCE
        https://pymolwiki.org/index.php/CGO_Shapes
    """
    from math import sin, cos, pi

    sides, quiet = int(sides), int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)

    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))
    rotmat = cpv.rotation_matrix(rotation, cpv.normalize(normal))

    def add_normal(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.NORMAL] + _xyz)

    def add_vertex(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.VERTEX] + cpv.add(center, _xyz))

    # Bottom
    obj.append(cgo.BEGIN)
    obj.append(cgo.TRIANGLE_FAN)
    if color:
        obj.extend([cgo.COLOR, *color])
    add_normal((0, 0, -1))
    add_vertex((0, 0, -height * 0.5))
    for i in range(sides + 1)[::-1]:
        angle = 2 * pi * i / sides
        x1, y1 = cos(angle) * radius, sin(angle) * radius
        add_vertex((x1, y1, -height * 0.5))
    obj.append(cgo.END)

    # Sides
    for i in range(sides + 1):
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLES)
        if color:
            obj.extend([cgo.COLOR, *color])
        angle1 = 2 * pi * i / sides
        x1, y1 = cos(angle1) * radius, sin(angle1) * radius
        angle2 = 2 * pi * (i + 1) / sides
        x2, y2 = cos(angle2) * radius, sin(angle2) * radius
        v1 = cpv.sub((x1, y1, -height * 0.5), (0, 0,  height * 0.5))
        v2 = cpv.sub((x2, y2, -height * 0.5), (0, 0,  height * 0.5))
        add_normal(cpv.cross_product(v1, v2))
        add_vertex((x1, y1, -height * 0.5))
        add_vertex((x2, y2, -height * 0.5))
        add_vertex((0, 0,  height * 0.5))
        obj.append(cgo.END)

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def torus(center=(0, 0, 0), normal=(0, 0, 1), radius=0.5, 
          cradius=.25, color="", samples=20, *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO torus object.
    REFERENCE
        cgobuilder, Copyright (c) Schrodinger, LLC
    """
    from math import sin, cos, pi

    samples, quiet = int(samples), int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)
    obj = []

    axis = cpv.cross_product(normal, (0., 0., 1.))
    angle = -cpv.get_angle(normal, (0., 0., 1.))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))

    def add_normal(xyz):
        return obj.extend(
            [cgo.NORMAL] + cpv.transform(matrix, xyz)
        )

    def add_vertex(xyz):
        return obj.extend(
            [cgo.VERTEX] + cpv.add(center, cpv.transform(matrix, xyz))
        )

    r = radius
    cr = cradius
    rr = 1.5 * cr
    dv = 2 * pi / samples
    dw = 2 * pi / samples
    v = 0.0
    w = 0.0

    while w < 2 * pi:
        v = 0.0
        c_w = cos(w)
        s_w = sin(w)
        c_wdw = cos(w + dw)
        s_wdw = sin(w + dw)

        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLE_STRIP)
        if color:
            obj.append(cgo.COLOR)
            obj.extend(color)

        while v < 2 * pi + dv:
            c_v = cos(v)
            s_v = sin(v)
            c_vdv = cos(v + dv)
            s_vdv = sin(v + dv)
            add_normal(
                [(r + rr * c_v) * c_w - (r + cr * c_v) * c_w,
                (r + rr * c_v) * s_w - (r + cr * c_v) * s_w,
                (rr * s_v - cr * s_v)]
            )
            add_vertex(
                [(r + cr * c_v) * c_w,
                (r + cr * c_v) * s_w,
                cr * s_v]
            )
            add_normal(
                [(r + rr * c_vdv) * c_wdw - (r + cr * c_vdv) * c_wdw,
                (r + rr * c_vdv) * s_wdw - (r + cr * c_vdv) * s_wdw,
                rr * s_vdv - cr * s_vdv]
            )
            add_vertex(
                [(r + cr * c_vdv) * c_wdw,
                (r + cr * c_vdv) * s_wdw,
                cr * s_vdv]
            )
            v += dv
        obj.append(cgo.END)
        w += dw

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj


@cmd.extend
def octahedron(center=(0, 0, 0), normal=(0, 0, 1), rotation=0,
               length=1.0, color="", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Return a CGO octahedron object.
    USAGE 
        octahedron [ center [, normal [, rotation, [, length [, color ]]]]]
    ARGUMENTS
        center = float3: object position  {default: 0, 0, 0}
        normal = float3: object orientation {default: 0, 0, 1}
        rotation = float: rotation around normal [rad] {default: 0}
        length = float: object side length {default: 1}
        color = string: object color {default: None}
    """
    from math import sqrt
    quiet = int(quiet)
    if _self.is_string(center):
        center = _self.safe_list_eval(center)
    if _self.is_string(normal):
        normal = _self.safe_list_eval(normal)
    if _self.is_string(length):
        length = _self.safe_list_eval(length)
    if color and isinstance(color, str):
        color = _self.get_color_tuple(color)

    obj = []

    axis = cpv.cross_product(normal, (0, 0, 1))
    angle = -cpv.get_angle(normal, (0, 0, 1))
    matrix = cpv.rotation_matrix(angle, cpv.normalize(axis))
    rotmat = cpv.rotation_matrix(rotation, cpv.normalize(normal))

    def add_normal(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.NORMAL] + _xyz)

    def add_vertex(xyz):
        _xyz = cpv.transform(rotmat, cpv.transform(matrix, xyz))
        return obj.extend([cgo.VERTEX] + cpv.add(center, _xyz))

    x, y, z = cpv.scale((0.5, 0.5, 1./sqrt(2.)), length)

    # Define vertices
    vertices = [
        [ x,  y,  0],  # 0
        [-x, -y,  0],  # 1
        [-x,  y,  0],  # 2
        [ x, -y,  0],  # 3
        [ 0,  0,  z],  # 4 (top)
        [ 0,  0, -z],  # 5 (bottom)
    ]

    # Define triangular faces
    faces = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [0, 5, 2],
        [2, 5, 1],
        [1, 5, 3],
        [3, 5, 0],
    ]

    for face in faces:
        v0, v1, v2 = [vertices[i] for i in face]
        r1 = cpv.sub(v1, v0)
        r2 = cpv.sub(v2, v0)
        obj.append(cgo.BEGIN)
        obj.append(cgo.TRIANGLES)
        if color:
            obj.extend([cgo.COLOR, *color])
        add_normal(cpv.cross_product(r1, r2))
        add_vertex(v0)
        add_vertex(v1)
        add_vertex(v2)
        obj.append(cgo.END)

    if not quiet:
        name = _self.get_unused_name("shape")
        _self.load_cgo(obj, name, zoom=0)

    return obj
