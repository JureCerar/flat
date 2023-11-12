from pymol import cmd

# TODO: Check for if it works for multi chain alignments
def alignment_mapping(seq1, seq2):
    """
    DESCRIPTION
        Returns an iterator with seq1 indices mapped to seq2 indices
    EXAMPLE
        >>> mapping = dict(alignment_mapping(s1, s2))
    """
    i, j = -1, -1
    for a, b in zip(seq1, seq2):
        if a not in ("-", "/"):
            i += 1
        if b not in ("-", "/"):
            j += 1
        if a != "-" and b != "-":
            yield i, j

def alignment_read(filename, format=None):
    """
    DESCRIPTION
        Wrapper for Bio.AlignIO.read that guesses alignment file format.
    """
    from Bio import AlignIO

    if not format:
        with open(filename) as handle:
            for line in handle:
                if len(line.rstrip()) > 0:
                    break
        if line.startswith("CLUSTAL") or line.startswith("MUSCLE"):
            format = "clustal"
        elif line.startswith(">P1;"):
            format = "pir"
        elif line.startswith(">"):
            format = "fasta"
        elif line.startswith("# STOCKHOLM"):
            format = "stockholm"
        elif line.startswith("Align "):
            format = "fatcat"
        elif line.startswith("# ProSMART Alignment File"):
            format = "prosmart"
        else:
            format = "emboss"

    with open(filename) as handle:
        return AlignIO.read(handle, format)




