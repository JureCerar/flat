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

from pymol import cmd, CmdException, cgo
import collections
import numpy as np


_hydro_moment = {
    # D. Bandyopadhyay, E.L. Mehler, Proteins, 2007, 72(2), 646–659. doi:10.1002/prot.21958 
    "bandyopadhyay-methler": {
        "CYS":  1.15, "ILE":  0.97, "LEU":  0.87, "PHE":  0.85, "VAL":  0.83,
        "TRP":  0.67, "TYR":  0.60, "MET":  0.54, "ALA":  0.33, "PRO":  0.32,
        "HIS":  0.25, "THR":  0.21, "SER":  0.05, "ARG": -0.01, "GLN": -0.05,
        "ASN": -0.07, "ASP": -0.22, "GLU": -0.24, "LYS": -0.40, "GLY":  0.00
    },
    # S.D. Black, D.R. Mould, Anal. Chem., 1991, 193, 72-82. doi:10.1016/0003-2697(91)90045-u
    "black-mould": {
        "ALA":  0.616, "ARG":  0.000, "ASN":  0.236, "ASP":  0.028, "CYS":  0.680,
        "GLN":  0.251, "GLU":  0.043, "GLY":  0.501, "HIS":  0.165, "ILE":  0.943,
        "LEU":  0.943, "LYS":  0.283, "MET":  0.738, "PHE":  1.000, "PRO":  0.711,
        "SER":  0.359, "THR":  0.450, "TRP":  0.878, "TYR":  0.880, "VAL":  0.825,
        "ASX":  0.132, "GLX":  0.147, "CYS2": 0.721
    },
    # D. Eisenberg, et. al, Faraday Symp. Chem. Soc., 1982, 17, 109-120. doi:10.1039/fs9821700109
    "eisenberg": {
        "ILE":  0.73, "PHE":  0.61, "VAL":  0.54, "LEU":  0.53, "TRP":  0.37,
        "MET":  0.26, "ALA":  0.25, "GLY":  0.16, "CYS":  0.04, "TYR":  0.02,
        "PRO": -0.07, "THR": -0.18, "SER": -0.26, "HIS": -0.40, "GLU": -0.62,
        "ASN": -0.64, "GLN": -0.69, "ASP": -0.72, "LYS": -1.10, "ARG": -1.80
    },
    # T. Jain, et. al, Bioinformatics, 2017, 33(23), 3758–3766. doi:10.1093/bioinformatics/btx519
    "jain": {
        "ALA":  0.064, "ARG": -0.358, "ASN":  0.143, "ASP": -0.476, "CYS":  0.610,
        "GLN":  0.516, "GLU": -0.807, "GLY":  0.000, "HIS":  0.028, "ILE":  1.724,
        "LEU":  1.720, "LYS": -1.198, "MET":  0.551, "PHE":  2.764, "PRO":  0.486,
        "SER":  0.076, "THR":  0.181, "TRP":  3.132, "TYR":  2.059, "VAL":  1.085
    },
    # J. Kyte, R.F. Doolittle, Journal of Molecular Biology, 1982, 157(1). 105–132, doi:10.1016/0022-2836(82)90515-0
    "kyte-doolittle": {
        "ILE":  4.5, "VAL":  4.2, "LEU":  3.8, "PHE":  2.8, "CYS":  2.5,
        "MET":  1.9, "ALA":  1.8, "GLY": -0.4, "THR": -0.7, "TRP": -0.9,
        "SER": -0.8, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
        "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
    },
    # J.L. Meek, PNAS, 1980, 77(3), 1632–1636. doi:10.1073/pnas.77.3.1632 
    "meek": {
        "TRP": 14.9, "PHE": 13.2, "ILE": 13.9, "LEU":  8.8, "TYR":  6.1,
        "MET":  4.8, "VAL":  2.7, "PRO":  6.1, "THR":  2.7, "ARG":  0.8,
        "ALA":  0.5, "GLY":  0.0, "HIS": -3.5, "CYS": -6.8, "LYS":  0.1,
        "SER":  1.2, "ASN":  0.8, "GLN": -4.8, "ASP": -8.2, "GLU": -16.9
    },
    # S. Miyazawa, R.L. Jernigan, Macromolecules, 1985, 18(3), 534–552. doi:10.1021/ma00145a039
    "miyazawa": {
        "ALA": 5.330, "ARG": 4.180, "ASN": 3.710, "ASP": 3.590, "CYS": 7.930,
        "GLN": 3.870, "GLU": 3.650, "GLY": 4.480, "HIS": 5.100, "ILE": 8.830,
        "LEU": 8.470, "LYS": 2.950, "MET": 8.950, "PHE": 9.030, "PRO": 3.870,
        "SER": 4.090, "THR": 4.490, "TRP": 7.660, "TYR": 5.890, "VAL": 7.630
    },
    # G. Rose, et. al, Science, 1985, 229(4716), 834–838. doi:10.1126/science.4023714 
    "rose": {
        "ALA": 0.74, "ARG": 0.64, "ASN": 0.63, "ASP": 0.62, "CYS": 0.91,
        "GLN": 0.62, "GLU": 0.62, "GLY": 0.72, "HIS": 0.78, "ILE": 0.88,
        "LEU": 0.85, "LYS": 0.52, "MET": 0.85, "PHE": 0.88, "PRO": 0.64,
        "SER": 0.66, "THR": 0.70, "TRP": 0.85, "TYR": 0.76, "VAL": 0.86,
    },
    # W.C. Wimley, S.H. White, Nature Structural & Molecular Biology, 1996, 3(10), 842–848. doi:10.1038/nsb1096-842 
    "wimley-white": {
        "ILE": -0.31, "LEU": -0.56, "PHE": -1.13, "VAL":  0.07, "MET": -0.23,
        "PRO":  0.45, "TRP": -1.85, "HIS":  0.17, "THR":  0.14, "GLUP": -0.01,
        "GLN":  0.58, "CYS": -0.24, "TYR": -0.94, "ALA":  0.17, "SER":   0.13,
        "ASN":  0.42, "ASPP": -0.07, "ARG":  0.81, "GLY":  0.01, "HISP":  0.96,
        "GLU":  2.02, "LYS":   0.99, "ASP":  1.23
    },
    # D. Eisenberg, A.D. McLachlan, Nature, 1986, 319(6050), 199–203. doi:10.1038/319199a0 
    "eisenberg-dg": {
        "GLY":  0.00, "ALA":  0.67, "VAL":  1.50, "LEU":  1.90, "ILE":  1.90,
        "PRO":  1.20, "CYS":  0.38, "MET":  2.40, "THR":  0.52, "SER":  0.01,
        "PHE":  2.30, "TRP":  2.60, "TYR":  1.60, "ASN": -0.60, "GLN": -0.22,
        "ASP": -1.20, "GLU": -0.67, "HIS":  0.64, "LYS": -0.57, "ARG": -2.10
    },
}


@cmd.extend
def hydropathy(selection="(all)", method="Black-Mould", vis=1, *, _self=cmd):
    """
    DESCRIPTION
        Assign hydrophobicity to amino-acid residues in selection.
    USAGE
        hydropathy [ selection [, method [, vis ]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        method = str: Hydrophobic scale. {default: 'Black-Mould'}
        vis = bool: Visualize results {default: True}
    """
    vis = int(vis)

    if not method.lower() in _hydro_moment:
        raise CmdException(f"Unknown method: '{method}'")
    
    hpi = _hydro_moment.get(method.lower())
    _self.alter(selection, "b = hpi.get(resn, 0.0)", space=locals())

    if vis:
        palette = ["tv_green", "white", "orange"]
        obj = _self.get_object_list(selection)[0]
        vmin, vmax = min(hpi.values()), max(hpi.values())
        _self.spectrum("b", " ".join(palette), selection, minimum=vmin, maximum=vmax)
        _self.ramp_new("hydrophobicity", obj, [vmin, vmax], palette)

    return


@cmd.extend
def hydropathy_moment(selection="(all)", method="Black-Mould", state=1, vis=1, *, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Get first-order hydrophobic moment vector for selection.
    USAGE
        hydropathy_moment [ selection [, method [, state [, vis ]]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        method = str: Hydrophobic scale. {default: 'Black-Mould'}
        vis = bool: Visualize result {default: True}
    LITERATURE
        B.D. Silverman, PNAS, 2001, 98(9), 4996-5001. doi:10.1073/pnas.081086198
    """
    vis, state, quiet = int(vis), int(state), int(quiet)

    if len(_self.get_object_list(selection)) != 1:
        raise CmdException("Multiple object in selection.")

    # Get list of residues
    atoms = _self.get_model(selection).atom
    residues = set(
        (a.segi, a.chain, a.resn, a.resi) for a in atoms)
    n_residues = len(residues)

    # Get dictionary of hydrophobic indices
    if not method.lower() in _hydro_moment:
        raise CmdException(f"Unknown method: '{method}'")
    hydro_dict = _hydro_moment.get(method.lower()).copy()

    # Assign hydrophobic index
    hydro = list()
    for segi, chain, resn, resi in residues:
        if resn not in hydro_dict.keys():
            print(
                f"Warning: Cannot assign hydrophobic index to residue: '{resn}'")
        # Default hydrophobicity value is zero?
        hydro.append(hydro_dict.get(resn, 0.0))
    hydro = np.array(hydro)

    # Get residues and molecule center of mass
    com = _self.centerofmass(selection, state=state)
    xyz = list()
    for segi, chain, resn, resi in residues:
        xyz.append(
            _self.centerofmass(f"//{segi}/{chain}/{resi}", state))
    xyz = np.array(xyz)

    # Get residue SASA
    tmp = _self.get_unused_name("_tmp")
    _self.create(tmp, selection, state, 1, zoom=0, quiet=1)
    try:
        _self.set("dot_solvent", 1, tmp)
        _self.set("dot_density", 4, tmp)
        sasa_dict = collections.defaultdict(float)
        _self.get_area(tmp, 1, load_b=1)
        _self.iterate(tmp, "sasa_dict[segi,chain,resi]+=b", space=locals())
        sasa = list()
        # Retain correct order of variables
        for segi, chain, resn, resi in residues:
            sasa.append(sasa_dict[segi, chain, resi])
        sasa = np.array(sasa)
    finally:
        _self.delete(tmp)

    # Calculate hydrophobic moment
    hydro_moment = np.zeros(3)
    for i in range(n_residues):
        hydro_moment += hydro[i] * sasa[i] * (xyz[i] - com)

    if not quiet:
        print(
            " Util: hydro_moment =",
            f"{np.linalg.norm(hydro_moment):.3f}",
            np.array2string(hydro_moment, precision=2)
        )

    if vis:
        color1 = _self.get_color_tuple("green")
        color2 = _self.get_color_tuple("orange")

        # Scale arrow shape
        radius = 0.5
        hlength = radius * 3.0
        hradius = hlength * 0.6

        vmin, vmax = _self.get_extent(selection, state)
        norm = np.linalg.norm(np.array(vmax) - np.array(vmin))
        v = hydro_moment / np.linalg.norm(hydro_moment) * norm / 2

        xyz1 = com - v
        xyz2 = com + v
        xyz3 = xyz2 + v / np.linalg.norm(v) * hlength

        name = _self.get_unused_name("hydro_moment")
        obj = [
            cgo.CYLINDER, *xyz1, *xyz2, radius, *color1, *color2,
            cgo.CONE, *xyz2, *xyz3, hradius, 0.0, *color2, *color2, 1.0, 0.0,
        ]
        _self.load_cgo(obj, name, state=state, zoom=1)

    return hydro_moment


# Autocomplete
cmd.auto_arg[0].update({
    "hydropathy": cmd.auto_arg[0]["zoom"],
    "hydropathy_moment": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "hydropathy": [cmd.Shortcut(_hydro_moment.keys()), "method", ""],
    "hydropathy_moment": [cmd.Shortcut(_hydro_moment.keys()), "method", ""],
})
