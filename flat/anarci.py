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

from pymol import cmd, CmdException
import tempfile
import subprocess
from pathlib import Path


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


def _parse_anarci(selection="all", scheme="imgt", exe="ANARCI", *, _self=cmd) -> dict:
    """Querry ANARCI for sequence in selection and return results""" 
    # Find and validate ANARCI executable
    if exe:
        exe = _self.exp_path(exe)
    else:
        import shutil
        exe = shutil.which("ANARCI")

    result = dict()

    with tempfile.TemporaryDirectory() as tmpdir:
        # Save selection sequence
        fasta = Path.joinpath(Path(tmpdir), "sequence.fasta")
        _self.save(fasta, selection)
        # Generate ANARCI report in CSV format
        output = Path.joinpath(Path(tmpdir), "out")
        try:
            subprocess.check_call([exe, "-i", str(fasta), "--scheme", scheme,
                                  "--csv", "-o", str(output),], shell=True)
        except subprocess.SubprocessError as e:
            raise CmdException(e)
        # Process csv reports
        for file in Path(tmpdir).glob("*.csv"):
            with open(file) as f:
                lines = f.read().splitlines()
            # Parse header and data
            header = lines[0].split(",")
            for line in lines[1:]:
                data = line.split(",")
                # Pack header infor into dictionary
                anarci_data = {k: v for k, v in zip(header[:13], data[:13])}
                # Parse antibody numbering  
                numbering = dict()
                for key, val in zip(header[13:], data[13:]):
                    if key[-1].isalpha():
                        resi, pos = int(key[:-1]), key[-1]
                    else:
                        resi, pos = int(key), ""
                    numbering[(resi, pos)] = val
                # Join data and save result
                anarci_data["numbering"] = numbering
                obj, chain = anarci_data.get("Id").split("_")
                result[(obj, chain)] = anarci_data

    return result


@cmd.extend
def anarci(selection="(all)", scheme="imgt", exe="ANARCI", *, quiet=0, _self=cmd):
    """
    DESCRIPTION
        Antibody Numbering and Antigen Receptor ClassIfication.
    USAGE
        anarci [ selection [, scheme [, exe ]]]
    ARGUMENTS
        selection = str: Atom selection. {default: all}
        scheme = str: Which numbering scheme should be used. {default: imgt}
        exe = str: Path to ANARCI executable. {default: }
    SCHEMES
        imgt, kabat, chothia, martin, aho, wolfguy
    SEE ALSO
        https://github.com/oxpig/ANARCI
    """
    # NOTE: It would be nice to call directly the ANARCI module, but on Windows
    # machines ANARCI does not work (no hmmer). I made it work via WSL wrapper
    # but I must then call an executable script and not module. It's a pain in the
    # ass but we work with what we have.

    # TODO:
    # - [ ] Implement residue renaming based on ANARCI output
    # - [ ] Figure our how to deal with really short sequences

    # Load CDR numbering scheme
    if scheme not in _CDR_DEFINITION.keys():
        raise CmdException(f"Unknown scheme: '{scheme}'")
    L1, L2, L3 = _CDR_DEFINITION[scheme]["L"]
    H1, H2, H3 = _CDR_DEFINITION[scheme]["H"]

    # Get ANARCI numbering for selection
    anarci_data = _parse_anarci(selection, scheme, exe)

    def within(num, lower, upper):
        """Check if number is in range"""
        return lower <= num <= upper

    # Loop over each chain in object
    for object in _self.get_names("public_nongroup_objects", 0, selection):
        # Create a selection group
        _self.delete(f"{object}_CDR")
        _self.group(f"{object}_CDR")
        for chain in _self.get_chains(object):
            # Get ANARCI numbering for selection
            df = anarci_data.get((object, chain), None)
            if not df:
                continue

            # Select CDR regions based on ANARCI output
            cdr_seq = ["", "", ""]
            chain_type = df.get("chain_type")
            numbering = df.get("numbering")
            if chain_type in ("H",):
                chain_name = "H" # Heavy chain 
                for resi, pos in numbering:
                    resn = numbering.get((resi, pos))
                    if resn == "-":
                        continue
                    if within(resi, *H1):
                        cdr_seq[0] += resn
                    elif within(resi, *H2):
                        cdr_seq[1] += resn
                    elif within(resi, *H3):
                        cdr_seq[2] += resn
            elif chain_type in ("K", "L",):
                chain_name = "L" # Light chain 
                for resi, pos in numbering:
                    resn = numbering.get((resi, pos))
                    if resn == "-":
                        continue
                    if within(resi, *L1):
                        cdr_seq[0] += resn
                    elif within(resi, *L2):
                        cdr_seq[1] += resn
                    elif within(resi, *L3):
                        cdr_seq[2] += resn

            # Print header information
            if not quiet:
                keys = ["hmm_species", "chain_type", "e-value", "score"]
                print(df['Id'])
                print("#", "|".join(keys))
                print("#","|".join(df[k] for k in keys))
                for i in range(3):
                    if not cdr_seq[i]:
                        continue
                    print(f"{chain_name}{i+1}:", cdr_seq[i])

            # Select and group CDR regions
            for i in range(3):
                if not cdr_seq[i]:
                    continue
                selection_name = f"{object}_{chain_name}{i+1}"
                _self.select(
                    selection_name,
                    f"{object} & c. {chain} & pepseq {cdr_seq[i]}",
                    merge=1,
                    enable=0,
                )
                _self.group(f"{object}_CDR", selection_name)
           
    return


# Autocomplete
cmd.auto_arg[0].update({
    "anarci": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "anarci": [cmd.Shortcut(_CDR_DEFINITION.keys()), "scheme", ""],
})
