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
import tempfile
import subprocess
import os

from . import one_letter


_CDR_DEFINITION = {
    "imgt": {
        # https://doi.org/10.1016/j.dci.2004.07.003
        "L": [(27, 38), (56, 65), (105, 117)],
        "H": [(27, 38), (56, 65), (105, 117)],
    },
    "kabat": {
        # https://doi.org/10.1084/jem.132.2.211 
        "L": [(24, 34), (50, 56), (89, 97)],
        "H": [(31, 35), (50, 65), (95, 102)],
    },
    "chothia": {
        # https://doi.org/10.1016/0022-2836(87)90412-8
        "L": [(26, 32), (50, 52), (91, 96)],
        "H": [(26, 32), (52, 56), (95, 102)],
    },
    "martin": {
        # https://doi.org/10.1016/j.molimm.2008.05.022
        "L": [(24, 34), (50, 56), (89, 97)],
        "H": [(26, 35), (50, 58), (95, 102)],
    },
    "aho": {
        # https://doi.org/10.1006/jmbi.2001.4662 
        "L": [(25, 40), (58, 77), (109, 137)],
        "H": [(25, 40), (58, 77), (109, 137)],
    },
    "wolfguy": {
        # https://doi.org/10.1002/prot.24756
        "L": [(551, 599), (651, 699), (751, 799)],
        "H": [(151, 199), (251, 299), (351, 399)],
    },
}


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
def anarci(selection="(all)", scheme="chothia", name="CDR", filename=None, exe="ANARCI.exe", *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Antibody Numbering and Antigen Receptor ClassIfication.
    USAGE
        anarci [ selection [, scheme [, name [, filename [, exe ]]]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        scheme = str: Which numbering scheme should be used. {default: chothia}
        name = str: a unique name for the selection. {default: CDR}
        filename = str: Save ANARCI output to file. {default: None}
        exe = str: Path to ANARCI executable. {default: }
    SCHEMES
        imgt, kabat, chothia, martin, aho, wolfguy
    SEE ALSO
        https://github.com/oxpig/ANARCI
    """
    # NOTE: It would be nice to call directly the ANARCI module, but on Windows
    # machines ANARCI does not work (no hmmer). I made it work via WSL wrapper
    # but I must then call an executable script and not module.

    # TODO: Implement residue renaming based on ANARCI output

    # TODO: Figure our how to deal with really short sequences

    # TODO: Implement saving ANARCI output to file

    # TODO: implement passing other variables to anarci
    # --restrict {ig,tr,heavy,light,H,K,L,A,B}
    # --assign_germline
    # --use_species {human,mouse,rat,rabbit,rhesus,pig,alpaca,cow}

    if not name:
        raise CmdException("Name cannot be empty")

    # Load CDR numbering scheme
    if scheme not in _CDR_DEFINITION.keys():
        raise CmdException(f"Unknown scheme: '{scheme}'")
    L1, L2, L3 = _CDR_DEFINITION[scheme]["L"]
    H1, H2, H3 = _CDR_DEFINITION[scheme]["H"]

    # Find and validate ANARCI executable
    if exe:
        exe = _self.exp_path(exe)
    else:
        import shutil
        exe = shutil.which("ANARCI")

    try:
        proc = subprocess.call([exe, "-h"], stdout=subprocess.DEVNULL,
                               stderr=subprocess.STDOUT, shell=True, text=True,)
        if proc < 0:
            raise CmdException(f"Broken executable: '{exe}'")
    except OSError as e:
        raise CmdException(f"Cannot execute '{exe!r}'") from e

    def irange(start, stop):
        """Inclusive range generator"""
        i = start
        while i <= stop:
            yield i
            i += 1
    
    # Remove named selection if it already exists
    _self.delete(f"{name}_H*")
    _self.delete(f"{name}_L*")

    # Loop over each chain in object
    for object in _self.get_names("public_nongroup_objects", 0, selection):
        for chain in _self.get_chains(object):

            # Get amino acid sequence
            residues = []
            _self.iterate(f"{object} & c. {chain} & guide",
                          "residues.append((resn, resv))", space=locals(),)
            seq = "".join(one_letter.get(x, "-") for x, _ in residues)

            # Process seqence with ANARCI
            proc = subprocess.run([exe, "--scheme", scheme, "-i", seq],
                                  shell=True, capture_output=True, text=True)
            if proc.returncode < 0:
                print(proc.stderr)
                raise CmdException(f"Process returned code: {proc.returncode}")

            # Parse ANARCI output
            lines = proc.stdout.splitlines()
            header = []
            cdr = ["", "", ""]
            for line in lines:
                # Header lines
                if line.startswith("#"):
                    header.append(line)
                    continue

                # End of entry
                if line.startswith("//"):
                    break

                # Process line
                col = line.split()
                if len(col) == 3:
                    region, resi, resn = col
                elif len(col) == 4:
                    region, resi, pos, resn = col

                # Skip missing residues
                if resn == "-":
                    continue
                resi = int(resi)

                # Determine CDR sequence
                if region == "L":
                    if resi in irange(*L1):
                        cdr[0] += resn
                    elif resi in irange(*L2):
                        cdr[1] += resn
                    elif resi in irange(*L3):
                        cdr[2] += resn
                elif region == "H":
                    if resi in irange(*H1):
                        cdr[0] += resn
                    elif resi in irange(*H2):
                        cdr[1] += resn
                    elif resi in irange(*H3):
                        cdr[2] += resn

            # Print header info
            if not quiet:
                print(15 * "-")
                print(f"|{object}|Chain {chain}|")
                for line in header:
                    print(line[1:].strip())

            # Select CDR regions
            for i in range(3):
                if not cdr[i]:
                    continue
                _self.select(f"{name}_{region}{i+1}",
                             f"{object} & c. {chain} & pepseq {cdr[i]}", merge=1, enable=0)

    return


# Autocomplete
cmd.auto_arg[0].update({
    "propka": cmd.auto_arg[0]["zoom"],
    "anarci": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "anarci": [cmd.Shortcut(_CDR_DEFINITION.keys()), "scheme", ""],
})
