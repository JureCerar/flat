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

from pymol import cmd
import numpy as np


@cmd.extend
def align_eig(selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with along it's eigen vector. 
    USAGE
        align_eig [ selection [, state ]]
    """
    xyz = np.array(_self.get_coords(selection, state)).T
    cov = np.cov(xyz)
    eval, evec = np.linalg.eig(cov)
    xyz = np.matmul(evec.T, xyz).T
    xyz -= np.mean(xyz, axis=0)
    _self.load_coords(xyz, selection, state=state)
    return


@cmd.extend
def align_3p(pk1="pk1", pk2="pk2", pk3="pk3", selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Align selection with along three selected points. 
    USAGE
        align_3p pk1, pk2, pk3, [ selection [, state ]]
    """
    v1 = np.array(*_self.get_coords(pk1, state))
    v2 = np.array(*_self.get_coords(pk2, state))
    v3 = np.array(*_self.get_coords(pk3, state))
    # Define a new base
    vx = v1 - v2
    vz = np.cross(v3 - v2, vx)
    vy = np.cross(vz, vx)
    # Normalize base
    vx /= np.linalg.norm(vx)
    vy /= np.linalg.norm(vy)
    vz /= np.linalg.norm(vz)
    # Stack to matrix
    base = np.array([vx, vy, vz]).reshape(3, 3)
    # Transform to new base
    # TODO: Account for center shift
    xyz = np.array(_self.get_coords(selection, state)).T
    xyz = np.matmul(base.T, xyz)
    _self.load_coords(xyz.T, selection, state=state)
    return
