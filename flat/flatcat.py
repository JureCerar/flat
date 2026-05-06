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

from pymol import cmd

@cmd.extend
def flatcat(*, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Checks if the cat is still flat
    USAGE
        flatcat 
    RETURNS
    : bool
        Returns `True` if the cat is flat
    """
    img = """
        |\---/| 
        | , , |      Yup, still flat...
        \_ x _/-..----.
     ___/  `    ,""-   \ 
    (__...'    __\      |'.___.':
      (_,...' (_,. __)./ '.....' 
    """
    if not quiet:
        print(img)
    return img