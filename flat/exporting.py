from pymol import cmd
import textwrap

one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

@cmd.extend
def get_pir(selection='(all)', state=-1, *, _self=cmd):
    # Default PIR settings
    chainbreak = "/"
    unknown = "-"
    width = 75

    # Get list of molecular objects
    objects = _self.get_names("public_objects", 0, selection)
    objects = [_ for _ in objects if _self.get_type(obj) == "object:molecule"]
    
    buffer = []
    for num, obj in enumerate(objects):
        # Extract info from model
        mdl = cmd.get_model(f"bca. {obj} & polymer")
        resn = [a.resn for a in mdl.atom]
        resi = [a.resi_number for a in mdl.atom]
        chain = [a.chain for a in mdl.atom]

        # Construct header and info lines
        buffer.append("")
        buffer.append(f">P{num + 1};{obj}")
        buffer.append(
            f"structureX:{obj}:{resi[0]}:{chain[0]}:{resi[-1]}:{chain[-1]}:::-1.00:-1.00"
        )

        # Get 1-letter amino acid sequence
        seq = ""
        prev = None
        for name, i in zip(resn, resi):
            if i != prev:
                if prev is not None and i != prev + 1:
                    seq += chainbreak
                else:
                    seq += one_letter.get(name, unknown)
                prev = i
        seq += "*"
        buffer += textwrap.wrap(seq, width)

    return "\n".join(buffer)


try:
    from pymol.exporting import savefunctions
    savefunctions.setdefault("pir", get_pir)
except ImportError:
    pass
