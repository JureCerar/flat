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
import tempfile
import os

@cmd.extend
def propka(selection="all", state=0, filename=None, vis=1, optargs=[], *, _self=cmd):
    """
    DESCRIPTION
        Predicts the pKa values of ionizable groups in proteins and
        protein-ligand complexes based on the 3D structure using PROPKA.
    USAGE
        propka [ selection [, state [, file [, vis [, optargs ]]]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        state = int: Object state (0 for all states). {default: 0}
        filename = str: Save PROPKA output to file. {default: None}
        vis = bool: Visualize residue pKa values. {default: True}
        optargs = list: Optional arguments to be passed to PROPKA. {default: []}
    """
    # See: https://propka.readthedocs.io/en/latest/api.html
    # See: https://tuilab.github.io/tui/update/2017/02/21/propka.html
    import propka.run as pk

    state, vis = int(state), int(vis)

    temp = os.path.join(tempfile.mkdtemp(), "temp.pdb")
    _self.save(temp, selection, state)

    try:
        pka = pk.single(temp, optargs, None, False)
    except:
        pass

    if vis:
        _self.show(
            "lines",
            f"((byres ({selection})) & (sc.|(n. CA|n. N&r. PRO)))"
        )
        _self.set("label_size", 12.0)
        _self.set("sphere_scale", 0.25, "(all)")
        for group in pka.conformations["AVR"].groups:
            if group.residue_type in ["ASP", "GLU", "HIS", "TYRF", "LYS", "ARG", "N+", "C-"]:
                chain_id = "" if group.atom.chain_id == "_" else group.atom.chain_id
                label = "'pKa={:.2f}'".format(group.pka_value)
                atom = "///{}/{}`{}/{}".format(
                    chain_id,
                    group.atom.res_name,
                    group.atom.res_num,
                    group.atom.name,
                )
                _self.label(atom, label)
                _self.show("spheres", atom)

    if filename:
        pka.write_pka(filename)

    return pka


@cmd.extend
def anarci(selection="all", *, _self=cmd):
    raise NotImplemented()

# Autocomplete
cmd.auto_arg[0].update({
    "propka": cmd.auto_arg[0]["zoom"],
    "anarci": cmd.auto_arg[0]["zoom"],
})