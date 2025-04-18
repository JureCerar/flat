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


def needle_alignment(s1, s2):
    '''
    DESCRIPTION
        Does a Needleman-Wunsch Alignment of sequence s1 and s2 and
        returns a Bio.Align.MultipleSeqAlignment object.
    SOURCE 
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    '''
    from Bio.Align import PairwiseAligner, substitution_matrices, MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    blosum62 = substitution_matrices.load("BLOSUM62")
    missing_codes = ''.join(set('JUO-.?').difference(blosum62.alphabet))
    blosum62 = blosum62.select(blosum62.alphabet + missing_codes)
    aligner = PairwiseAligner(
        internal_open_gap_score=-10,
        extend_gap_score=-0.5,
        substitution_matrix=blosum62
    )
    alns = aligner.align(s1, s2)
    return MultipleSeqAlignment([
        SeqRecord(Seq(alns[0]), "s1"),
        SeqRecord(Seq(alns[1]), "s2"),
    ])
