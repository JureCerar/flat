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

import tempfile
import subprocess
from pathlib import Path
from pymol import cmd

from .exporting import save_xyzr

@cmd.extend
def msms(selection="polymer", state=1, density=3, name=None,
         filename=None, *, exe="msms", _self=cmd):
    """
    DESCRIPTION
        Run MSMS on the given selection and load the generated surface as a CGO.
    USAGE
        msms [ selection [, state [, density [, name [, filename ]]]]]
    ARGUMENTS
        selection = str: atom selection {default: polymer}
        state = int: object state {default: 1}
        density = float: MSMS surface point density. {default: 3}
        name = str: name of CGO object to create. {default: None}
        filename = str: Name of output file (without extension). {default: None}
    SOURCE
        From PSICO (c) 2010-2012 Thomas Holder, MPI for Developmental Biology
    """
    density = float(density)
    hdensity = density + 2
    solvent_radius = _self.get("solvent_radius")

    if exe:
        exe = _self.exp_path(exe)
    else:
        import shutil
        exe = shutil.which("msms")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_if = Path().joinpath(tmpdir, "tmp.xyzr")
        if not filename:
            filename = Path().joinpath(tmpdir, "tmp")
        save_xyzr(str(tmp_if), selection, state, _self=_self)
        subprocess.check_call([exe, "-density", str(density), "-hdensity", str(hdensity),
                               "-if", str(tmp_if), "-of", str(filename), "-no_area",
                               "-probe_radius", solvent_radius,], cwd=tmpdir, shell=True)
        if not name:
            name = _self.get_unused_name("surface")
        _self.load(str(filename) + ".face", name, _self=_self)


# Autocomplete
cmd.auto_arg[0].update({
    "msms": cmd.auto_arg[0]["zoom"],
})
