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

from setuptools import find_packages, setup

try:
  from flat import __version__
except:
  print("Could not import package version")
  __version__ = None

setup(
  name="flat",
  version=__version__,
  description="PyMol Script Library/Collection/Repository",
  author="Jure Cerar",
  license="GNU GPL 3.0",
  url="https://github.com/JureCerar/flat",
  packages=find_packages(),
  include_package_data=True,
  requires=["pymol"],
)
