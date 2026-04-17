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

"""
:mod:`flat.propka`
==================
Module for computing pKa values of ionizable groups using `PROPKA`_.

Installation
------------
You can easily install PROPKA with:

.. code:: bash

    pip install propka

.. PROPKA:
    https://propka.readthedocs.io/en/latest/
"""

from pymol import cmd

# TODO:
# - Add propka_plot function 

@cmd.extend
def propka(selection="all", state=0, filename=None, vis=1, optargs=[], *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Predicts the pKa values of ionizable groups in proteins and
        protein-ligand complexes based on the 3D structure using PROPKA.
    USAGE
        propka [ selection [, state [, file [, vis [, optargs ]]]]]
    ARGUMENTS
        selection : str, optional
            Atom selection.
        state : int, default = 0
            Object state (0 for all states).
        filename : str, optional
            Save PROPKA output to file.
        vis : bool, default = True
            Visualize residue pKa values.
        optargs : List[str], optional
            Optional arguments to be passed to PROPKA.
    RETURNS
        : propka.MolecularContainer
            Propka model for selection.
    REFERENCE
        https://propka.readthedocs.io/en/latest/
    """
    # See: https://propka.readthedocs.io/en/latest/api.html
    # See: https://tuilab.github.io/tui/update/2017/02/21/propka.html
    import propka.run as pk
    import tempfile
 
    # Save PDB file and process it with PROPKA
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        _self.save(tmp.name, selection, int(state))
        pka = pk.single(tmp.name, optargs, None, False)

    # Show ionizable groups
    if int(vis):
        _self.show(
            "lines",
            f"((byres ({selection})) & (sc.|(n. CA|n. N&r. PRO)))"
        )
        _self.set("label_size", 12.0)
        _self.set("sphere_scale", 0.25, "(all)")
        for group in pka.conformations["AVR"].groups:
            if group.residue_type in ["ASP", "GLU", "HIS", "TYR", "LYS", "ARG", "N+", "C-"]:
                chain_id = "" if group.atom.chain_id == "_" else group.atom.chain_id
                label = "'pKa={:.2f} ({:.1f})'".format(group.pka_value, group.model_pka)
                atom = "///{}/{}`{}/{}".format(
                    chain_id,
                    group.atom.res_name,
                    group.atom.res_num,
                    group.atom.name,
                )
                _self.label(atom, label)
                _self.show("spheres", atom)
    
    if not int(quiet):
        pi_folded, pi_unfolded = pka.get_pi()
        print(f"Util: pI = {pi_folded:.2f} ({pi_unfolded:.2f})")

    # Save file
    if filename:
        pka.write_pka(filename)

    return pka


# Autocomplete
cmd.auto_arg[0].update({
    "propka": cmd.auto_arg[0]["zoom"],
})
