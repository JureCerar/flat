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

from pymol import cmd, CmdException
import textwrap

@cmd.extend
def save_csv(filename, selection="(all)", var="b", mode="atom", quiet=1, *, _self=cmd):
    """
    DESCRIPTION
        Save property from selection to CSV file.
    USAGE
        save_csv file [, selection [, var [, mode ]]]
    ARGUMENTS
        file = str: Output file name.
        selection = str: Atom selection. {default: all}
        var = str: Property to save. {default: b}
        mode = str: Iterator mode: atom or res. {default: "atom"}
    """
    if len(_self.get_object_list(selection)) > 1:
        raise CmdException("Multiple objects in selection.")

    data = []
    if mode.lower() == "atom":
        # Iterate by atoms
        data.append(["chain", "resn", "resi", "name", "value"])
        _self.iterate(
            selection,
            "data.append([chain, resn, resi, name, {}])".format(var),
            space=locals(),
        )

    elif mode.lower() == "res":
        # Iterate by residues
        data.append(["chain", "resn", "resi", "value"])
        _self.iterate(
            f"bca. ({selection})",
            "data.append([chain, resn, resi, {}])".format(var),
            space=locals(),
        )

    else:
        raise CmdException(f"Invalid mode: '{mode}'")

    with open(filename, "w") as f:
        for line in data:
            print(*line, sep=",", file=f)

    if not int(quiet):
        print(f'Save: wrote "{filename}"')

    return


@cmd.extend
def save_ndx(filename, quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Save all selections as GROMACS index (.ndx) file.
    USAGE
        save_ndx filename [, quiet ]]
    """
    quiet = int(quiet)

    selections = cmd.get_names("public_selections")
    groups = dict()

    for name in selections:
        groups[name] = list()
        cmd.iterate(name, "group.append(index)", space={"group": groups[name]})

    with open(filename, "w") as handle:
        nlines = 15  # Number of indeces per line
        for name, group in groups.items():
            if not quiet:
                print(f"Saving group: '{name}' len(group atoms)")
            print(f"[ {name} ]", file=handle)
            for i in range(0, len(group), nlines):
                print(
                    " ".join(f"{x:5}" for x in group[i:i+nlines]), file=handle)

    return


@cmd.extend
def get_pir(selection='(all)', state=-1, *, _self=cmd):
    """
    DESCRIPTION
        Get sequence or sequence alignment in PIR format.
    USAGE
        get_pir [selection]
    """
    from pymol.exporting import _resn_to_aa as one_letter
    
    # Get list of molecular objects
    objects = _self.get_names("public_objects", 0, selection)
    objects = [_ for _ in objects if _self.get_type(obj) == "object:molecule"]
    
    buffer = []
    for obj in objects:
        # Extract info from model
        mdl = cmd.get_model(f"bca. {obj} & polymer")
        resn = [a.resn for a in mdl.atom]
        resi = [a.resi_number for a in mdl.atom]
        chain = [a.chain for a in mdl.atom]

        # Construct header and info lines
        buffer.append("")
        buffer.append(f">P1;{obj}")
        buffer.append(
            f"structureX:{obj}:{resi[0]}:{chain[0]}:{resi[-1]}:{chain[-1]}:::-1.00:-1.00"
        )

        # Get 1-letter amino acid sequence
        seq = ""
        prev = None
        for name, i in zip(resn, resi):
            if i != prev:
                if prev is not None and i != prev + 1:
                    seq += "/" # Chain break
                else:
                    seq += one_letter.get(name, "-")
                prev = i
        seq += "*"
        buffer += textwrap.wrap(seq, 75)

    return "\n".join(buffer)


try:
    from pymol.exporting import savefunctions
    savefunctions.setdefault("csv", save_csv)
    savefunctions.setdefault("ndx", save_ndx)
    savefunctions.setdefault("pir", get_pir)
except ImportError:
    pass
