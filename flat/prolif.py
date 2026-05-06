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
:mod:`flat.prolif`
==================
Module for analyzing interactions with `ProLIF`_.

Installation
------------
You can easily install ProLIF (and rdkit) with:

.. code:: bash

    pip install prolif rdkit

.. ProLIF:
    https://prolif.readthedocs.io
"""

from pymol import cmd

from .querying import iterate_state_to_list

# Custom colors for contacts
CONTACT_COLORS = {
    "VdWContact": "magenta",
    "Hydrophobic": "orange",
    "PiStacking": "gray",
    "HBDonor": "yellow",
    "HBAcceptor": "yellow",
    "Anionic": "red",
    "Cationic": "red",
    "CationPi": "purple",
    "PiCation": "purple",
}

@cmd.extend
def prolif(sele_lig, sele_pro="polymer", state=-1, prefix="prolif.", *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Find interactions with ProLIF.
    USAGE
        prolif sele_lig [, sele_pro [, state [, prefix ]]]
    ARGUMENTS
        sele_lig : str
            Ligand selection.
        sele_pro : str, default = 'polymer'
            Protein selection ()
        state : int, default = -1
            Object state (-1 for current state).
        prefix : str, default = 'prolif.'
            Name of output group.
    EXAMPLE
        >>> fetch 1eve
        >>> h_add
        >>> prolif resn E20
    SOURCE
        From PSICO (c) 2025 Thomas Holder
    SEE ALSO
        https://prolif.readthedocs.io
    """
    import prolif as plf
    from rdkit import Chem

    idx_pro = iterate_state_to_list(state, sele_pro,
                                    "f'{model}`{index}'", _self=_self)
    idx_lig = iterate_state_to_list(state, sele_lig,
                                    "f'{model}`{index}'", _self=_self)

    rdkit_mol_pro = Chem.MolFromPDBBlock(_self.get_str("pdb", sele_pro, state),
                                         removeHs=False, sanitize=False)
    rdkit_mol_lig = Chem.MolFromMolBlock(_self.get_str("mol", sele_lig, state),
                                         removeHs=False, sanitize=False)

    for rdkit_mol in [rdkit_mol_pro, rdkit_mol_lig]:
        failed_flag = Chem.SanitizeMol(rdkit_mol, catchErrors=True)
        if failed_flag:
            print(f"Sanitize failed for {rdkit_mol} on flag {failed_flag}")

    fp = plf.Fingerprint()
    fp.run_from_iterable([plf.Molecule(rdkit_mol_lig)],
                         plf.Molecule(rdkit_mol_pro), n_jobs=1)

    assert len(fp.ifp) == 1

    i_pro_all = set()

    if prefix.endswith("."):
        _self.group(prefix.removesuffix("."))

    pseudo = _self.get_unused_name("_pseudo")

    for inter_num, interactions in enumerate(fp.ifp[0].values(), 1):
        for inter_type, interaction_tuple in interactions.items():
            for inter_obj in interaction_tuple:
                indices = inter_obj["parent_indices"]

                if not int(quiet):
                    print(f"{inter_type:11s} {inter_obj['distance']:.2f} {indices}")

                i_lig = indices["ligand"]
                i_pro = indices["protein"]
                s_pro = " ".join(idx_pro[i] for i in i_pro)
                s_lig = " ".join(idx_lig[i] for i in i_lig)

                i_pro_all.update(i_pro)

                if len(i_lig) > 1:
                    _self.pseudoatom(pseudo, s_lig, resi=inter_num,
                                     name="LIG", state=state)
                    s_lig = f"/{pseudo}///`{inter_num}/LIG"

                if len(i_pro) > 1:
                    _self.pseudoatom(pseudo, s_pro, resi=inter_num,
                                     name="PRO", state=state)
                    s_pro = f"/{pseudo}///`{inter_num}/PRO"
                    
                objname = f"{prefix}{inter_type}"
                _self.distance(objname, s_pro, s_lig)

                color = CONTACT_COLORS.get(inter_type)
                if color:
                    _self.color(color, objname)
    
    _self.delete(pseudo)


# Autocomplete
cmd.auto_arg[0].update({
    "prolif": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "prolif": cmd.auto_arg[0]["zoom"],
})
