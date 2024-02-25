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

if __name__.endswith(".fullinit"):
    from . import init
   
    init()

else:
    from pathlib import Path
    from importlib import util
    import sys

    try:
        _script_path = __script__  # type: ignore[name-defined]
    except NameError:
        raise ImportError("invalid invocation of `flat.fullinit`") from None

    path = Path(_script_path).parent / "__init__.py"
    spec = util.spec_from_file_location("flat", path)
    assert spec is not None
    assert spec.loader is not None
    module = util.module_from_spec(spec)
    assert module is not None
    sys.modules["flat"] = module
    spec.loader.exec_module(module)

    import flat.fullinit
