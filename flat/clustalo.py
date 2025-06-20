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

from pymol import cmd
import tempfile
import subprocess 

# Formats supported by CLUSTAL 
_CLUSTAL_FMT = ["a2m", "fasta", "clustal", "msf", "phylip", "selex", "stockholm", "vienna"]


@cmd.extend
def clustalo(selection="all", fmt="clustal", exe="clustalo", *, _self=cmd):
    """
    DESCRIPTION
        Multiple sequence alignment with Clustal O
    USAGE
        clustalo [ selection [, fmt [, exe ]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        fmt = str: Output format. {default: clustal}
        exe = str: Path to CLUSTALO executable. {default: }
    """
    # TODO:
    # - [ ] Add result visualization

    # Find and validate CLUSTALO executable
    if exe:
        exe = _self.exp_path(exe)
    else:
        import shutil
        exe = shutil.which("clustalo")

    # Save sequence and process it with CLUSTALO
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        _self.save(tmp.name, selection)
        proc = subprocess.run([exe, "--outfmt", fmt, "-i", tmp.name],
                              shell=True, capture_output=True, text=True)
        
    for line in proc.stdout.splitlines():
        print(line)

    return


# Autocomplete
cmd.auto_arg[0].update({
    "clustalo": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "clustalo": [cmd.Shortcut(_CLUSTAL_FMT), "fmt", ""],
})
