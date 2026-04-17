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

"""
:mod:`flat.exporting`
=====================
Module that allows to save files in additional formats.
"""

from pymol import cmd
import textwrap
import csv


@cmd.extend
def save_gro(filename, selection="all", state=0, *, _self=cmd):
    """
    DESCRIPTION
        Save structure in GROMACS configuration (.gro) format.
    USAGE
        save_gro filename [, selection [, state ]]
    ARGUMENTS
        filename : str
            ... 
        selection : str, optional
            Atom selection.
        state : int, default = 0
            ...
    """
    # NOTE to self: f.write does some strange things?! Just use print.
    with open(filename, "w") as f:
        print("Generated with PyMOL", file=f)
        print(_self.count_atoms(selection, state=state), file=f)
        fmt = "{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}"
        callback = lambda *args: print(fmt.format(*args), file=f)
        _self.iterate_state(
            state,
            selection,
            "callback(int(resi), resn, name, index, x/10, y/10, z/10)",
            space=locals()
        )
        try:
            box = [x/10 for x in _self.get_symmetry(selection, state)[:3]]
        except:
            box = [0, 0, 0]
        print("{:10.5f}{:10.5f}{:10.5f}".format(*box), file=f)


@cmd.extend
def save_csv(filename, selection="all", selector=None, var="b", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Save property from selection to CSV file.
    USAGE
        save_csv filename [, selection [, selector [, var ]]]
    ARGUMENTS
        filename : str
            Output file name.
        selection : str, default = 'all'
            Atom selection.
        selector : List(str), default = None
            List or string of atom selectors to export properties. It
            defaults to: `[chain, resn, resi, name]`.
        var : str, default = 'b'
            Property to export.
    """
    if isinstance(selector, str):
        selector = selector.split()
    if not selector:
        selector = ["chain", "resn", "resi", "name"]
    selector_key = ",".join(selector)

    # NOTE: Load to dict to make sure output has one entry per selector; If there
    # are multiple atoms that fit selector just store last one and say it's a feature
    data = dict()
    _self.iterate(selection, f"data[{selector_key}]={var}", space=locals())
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f,)
        writer.writerow(selector + [var])
        for key, val in data.items():
            writer.writerow([*key, val])

    if not int(quiet):
        print(f"Save: Wrote {filename!r}")


@cmd.extend
def save_ndx(filename, *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Save all selections as GROMACS index (.ndx) format.
    USAGE
        save_ndx filename
    ARGUMENTS
        filename : str
            Output file name.
    """
    selections = _self.get_names("public_selections")
    groups = dict()

    # Gather all all indices
    for sele in selections:
        groups[sele] = list()
        _self.iterate(sele, "group.append(index)",
                      space={"group": groups[sele]})
        if not int(quiet):
            print(f"Saving group: {sele!r} ({len(groups[sele])} atoms)")

    NUM = 15  # Number of indices per line
    with open(filename, "w") as f:
        for name, group in groups.items():
            print(f"[ {name} ]", file=f)
            for i in range(0, len(group), NUM):
                print(" ".join(f"{x:5}" for x in group[i:i+NUM]), file=f)


@cmd.extend
def save_pir(filename, selection="all", *, _self=cmd):
    """
    DESCRIPTION
        Save sequence in Protein Information Resource (PIR) format.
    USAGE
        save_ndx filename [, selection ]
    ARGUMENTS
        filename : str
            Output file name.
        selection : str, default = 'all'
            Atom selection.
    """
    from pymol.exporting import _resn_to_aa as one_letter
    
    # Get list of molecular objects
    objects = _self.get_names("public_objects", 0, selection)
    objects = [o for o in objects if _self.get_type(o) == "object:molecule"]
    
    buffer = []
    for obj in objects:
        # Extract info from model
        mdl = _self.get_model(f"bca. {obj} & polymer")
        resn = [a.resn for a in mdl.atom]
        resi = [a.resi_number for a in mdl.atom]
        chain = [a.chain for a in mdl.atom]
        # Construct header and info lines
        buffer.extend(
            [
                "",
                f">P1;{obj}",
                f"structureX:{obj}:{resi[0]}:{chain[0]}:{resi[-1]}:{chain[-1]}:::-1.00:-1.00"
            ]
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

    with open(filename, "w") as f:
        print("\n".join(buffer), file=f)


@cmd.extend
def save_mda(filename, selection="all", *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Save a trajectory file with the MDAnalysis library.
    USAGE
        save_mda filename [, selection ]
    ARGUMENTS
        filename : str
            Output file name.
        selection : str, default = 'all'
            Atom selection.
    """
    import MDAnalysis
    import tempfile

    # Generate and read temporary topology
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
        _self.save(tmp, selection)
        tmp.flush()
        u = MDAnalysis.Universe(tmp, in_memory=True, trajectory=True)

    # if os.path.exists(filename):
    #     os.remove(filename)

    n_frames = _self.count_states(selection)
    with MDAnalysis.Writer(filename, u.atoms.n_atoms) as handle:
        for state in range(1, n_frames + 1):
            u.atoms.positions = _self.get_coords(selection, state)
            u.trajectory.dimensions = _self.get_symmetry(selection, state)[:6]
            handle.write(u)

    if not int(quiet):
        print(f"Wrote {n_frames} frames to file: {filename!r}")


@cmd.extend
def save_colored_fasta(filename, selection="all", invert=False, *, _self=cmd):
    """
    DESCRIPTION
        Save a HTML file with colored (by C-alpha atoms) FASTA sequence.
    USAGE
        save_colored_fasta filename [, selection [, invert ]]]
    ARGUMENTS
        filename : str
            Output file name.
        selection : str, default = 'all'
            Atom selection.
        invert : bool, default = False
            Color background instead of character.
    SOURCE
        From PSICO (c) 2012 Thomas Holder, MPI for Developmental Biology
    """
    # See: https://github.com/speleo3/pymol-psico/blob/master/psico/fasta.py
    from pymol import Scratch_Storage
    from pymol.exporting import _resn_to_aa as one_letter

    invert = int(invert)
    selection = f"({selection}) and guide"

    html = []
    stored = Scratch_Storage()

    def callback(resv, resn, color):
        """ Translate to residue to HTML."""
        if stored.resv is None:
            stored.resv = resv - (resv % 70)
        while stored.resv+1 < resv:
            callback(stored.resv+1, None, 222)
        stored.resv += 1
        if stored.resv % 70 == 1:
            html.append("<br>")

        rgb = cmd.get_color_tuple(color)
        rgb = [int(x*255) for x in rgb]
        aa = one_letter.get(resn, "?") if resn else "-"

        if not resn:
            color = "#000000"
            background = ""
        elif invert:
            cs = round((rgb[0] * 299 + rgb[1] * 587 + rgb[2] * 114) / 1000)
            color = "black" if cs > 125 else "white"
            background = "background-color:#{:02x}{:02x}{:02x}".format(*rgb)
        else:
            color = "#{:02x}{:02x}{:02x}".format(*rgb)
            background = ""

        html.append(
            f"<span style=color:{color};{background}>{aa}</span>"
        )
        if stored.resv % 10 == 0:
            html.append("<span> </span>")
        return

    for obj in _self.get_object_list(selection):
        for chain in _self.get_chains(f"model {obj} and ({selection})"):
            sele = f"m. {obj} & c. '{chain}' and ({selection})"
            html.append(f"&gt;{obj}|{chain}")
            stored.resv = 0
            _self.iterate(sele, "callback(resv, resn, color)", space=locals())
            html.append("<br>")

    with open(filename, "w") as handle:
        print("<html><body style='font-family:consolas'>", file=handle)
        print("".join(html), file=handle)
        print("</body></html>", file=handle)


@cmd.extend
def save_xyzr(filename, selection="all", state=1, *, _self=cmd):
    """
    DESCRIPTION
        Write the given selection to an xyzr or xyzrn (determined by extension)
        file for MSMS.
    USAGE
        save_xyzr filename [, selection [, state ]]]
    ARGUMENTS
        filename : str
            Output file name.
        selection : str, default = 'all'
            Atom selection.
        state : int, default = 1
            Object state.
    """
    if filename.endswith("xyzrn"):
        expr = "callback(x, y, z, vdw, name, resn, resi)"
        fmt = "%.3f %.3f %.3f %.2f 1 %s_%s_%s"
    else:
        expr = "callback(x, y, z, vdw)"
        fmt = "%.3f %.3f %.3f %.2f"
    with open(filename, "w") as f:
        callback = lambda *args: print(fmt % args, file=f)
        _self.iterate_state(state, selection, expr, space=locals())


# Register extensions
try:
    from pymol.exporting import savefunctions
    savefunctions.setdefault("gro", save_gro)
    savefunctions.setdefault("csv", save_csv)
    savefunctions.setdefault("ndx", save_ndx)
    savefunctions.setdefault("pir", save_pir)
    for ext in ["dcd", "xtc", "trj", "crd"]:
        savefunctions.setdefault(ext, save_mda)
    for ext in ["xyzr", "xyzrn",]:
        savefunctions.setdefault(ext, save_xyzr)
except ImportError:
    pass

# Autocomplete
cmd.auto_arg[1].update({
    "save_gro": cmd.auto_arg[1]["save"],
    "save_csv": cmd.auto_arg[1]["save"],
    "save_ndx": cmd.auto_arg[1]["save"],
    "save_pir": cmd.auto_arg[1]["save"],
    "save_mda": cmd.auto_arg[1]["save"],
    "save_colored_fasta": cmd.auto_arg[1]["save"],
    "save_xyzr": cmd.auto_arg[1]["save"],
})
