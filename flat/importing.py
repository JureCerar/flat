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
import collections
import re
import io
import csv

@cmd.extend
def load_all(pattern, group='', quiet=1, *, _self=cmd, **kwargs):
    '''
    DESCRIPTION
        Load all files matching given globbing pattern
    '''
    # import glob
    # import os
    # TODO: Idea structure first, everything else later
    # filenames = glob.glob(cmd.exp_path(pattern))
    # for filename in filenames:
    #     if not quiet:
    #         print(' Loading ' + filename)
    #     _self.load(filename, **kwargs)
    raise NotImplemented


@cmd.extend
def load_csv(filename, selection="(all)", prop=None, var="b", vis=0, *, _self=cmd):
    """
    DESCRIPTION
        Load property from CSV file to selection. Atom selector is automaticaly
        generated from CSV header.
    USAGE
        load_csv filename [, selection [, prop [, var [, vis ]]]]
    ARGUMENTS
        filename = str: Input file name.
        selection = str: Atom selection. {default: all}
        prop = str: Property name to import (first by default). {default: None}
        var = str: Variable where property is saved. {default: b}
        vis = int: Visualize output. {default: 1}
    """
    class FactoryCounter:
        def __init__(self, factory):
            self.factory = factory
            self.count = 0

        def __call__(self):
            self.count += 1
            return self.factory()

    with open(filename, "r") as f:
        lines = f.read()

    if not csv.Sniffer().has_header(lines):
        raise CmdException("CSV file does not contain header")
    
    reader = csv.DictReader(
        io.StringIO(lines),
        dialect=csv.Sniffer().sniff(lines),
    )

    # Valid PyMol selectors
    SELECTORS = ["object", "symbol", "name", "resn", "resi", "alt",
                 "chain", "segi", "flag", "id", "index", "ss",]
    selector = list()
    for fieldname in reader.fieldnames:
        if fieldname in SELECTORS:
            selector.append(fieldname)
    if not selector:
        raise ValueError("No selector rows found in CSV file")
    selector_key = ",".join(selector)

    if prop:
        if prop not in reader.fieldnames:
            raise ValueError("Property not found in CSV file")
    else:
        # Select first valid property
        for fieldname in reader.fieldnames:
            if fieldname not in SELECTORS:
                prop = fieldname
                break
        else:
            raise ValueError("No property rows found in CSV file")

    data = collections.defaultdict(factory := FactoryCounter(float))
    for row in reader:
        key = tuple(row[_] for _ in selector)
        data[key] = row[prop]

    _self.alter(selection, f"{var}=data[{selector_key}]", space=locals())

    # Give dictionary miss count
    if factory.count:
        print(f"Warning: {factory.count} atom(s) were not found in CSV file")

    if vis:
        palette = ["marine", "silver", "red"]
        obj = _self.get_object_list(selection)[0]
        values = list()
        _self.iterate(selection, f"values.append({var})", space=locals())
        vmin, vmax = min(values), max(values)
        _self.spectrum(var, " ".join(palette), selection, vmin, vmax)
        _self.ramp_new(f"{prop}_ramp", obj, (vmin, vmax), palette)

    return data


@cmd.extend
def load_topol(filename, selection="(all)", *, _self=cmd):
    """
    DESCRIPTION
        Load GROMACS topology (.top) file to selection.
    USAGE
        load_topol filename [, selection ]
    """
    atom_list = collections.defaultdict(lambda: [0, 0])
    current = None

    def get_keyword(string):
        """ Scan a string for keywords i.e. '[ keyword ]' """
        k = re.search(r'\[(.*?)\]', string)
        return k.string.split()[1] if k else None

    for line in open(filename, "r"):
        line = line.split(";")[0].strip()
        if not line:
            continue

        key = get_keyword(line)
        if key:
            current = key
            continue

        if current == "atoms":
            # nr, type, resnr, residue, atom, cgnr, pcharge, mass
            col = line.split()
            try:
                index, pcharge, mass = int(
                    col[0]), float(col[6]), float(col[7])
                atom_list[index] = [pcharge, mass]
            except:
                print(f"Cannot parse line: '{line}'")

    _self.alter(
        selection,
        "partial_charge, mass = atom_list[index]",
        space=locals()
    )
    return


@cmd.extend
def load_ndx(filename, quiet=0, *, _self=cmd):
    """
    DESCRIPTION
        Read a GROMACS index (.ndx) file as a selection.
    USAGE
        read_ndx filename [, quiet ]
    """
    quiet = int(quiet)
    groups = dict()

    def get_keyword(string):
        """ Scan a string for keywords i.e. '[ keyword ]' """
        k = re.search(r'\[(.*?)\]', string)
        return k.string.split()[1] if k else None

    current = None
    for line in open(filename, "r"):
        line = line.split(";")[0]
        key = get_keyword(line)
        if key:
            current = key
            groups[current] = []
        else:
            groups[current] += [int(i) for i in line.split()]

    for name, group in groups.items():
        # Pymol does not like '&' and '|' characters and reserved keywords
        name = name.replace("&", "_and_")
        name = name.replace("|", "_or_")
        if name.lower() in ["sidechain", "backbone"]:
            name = name + "_"

        if len(group) == 0:
            print(" Empty group: '%s'" % name)
            continue

        if not quiet:
            print("Loading group: '{}' ({} atoms)".format(name, len(group)))

        # Pymol does not like long selection strings so we select the group iteratively
        every = 15
        _self.select(f"({name})", f"index {group[0]}")
        for i in range(1, len(group), every):
            buffer = "index " + "+".join(str(x) for x in group[i:i+every])
            _self.select(name, f"({name}) | {buffer}", quiet=1)

    return


@cmd.extend
def load_smi(filename, oname="", discrete=-1, quiet=1, multiplex=None, zoom=-1, _self=cmd):
    """
    DESCRIPTION
        Load a SMILES file with an openbabel backend
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    """
    import tempfile
    import subprocess
    import os

    if not oname:
        oname = os.path.basename(filename).rsplit(".", 1)[0]

    outfile = tempfile.mktemp(".sdf")

    try:
        subprocess.check_call(["obabel", filename, "-O", outfile, "--gen3D"])
        _self.load(outfile, oname, discrete=discrete, quiet=quiet,
                   multiplex=multiplex, zoom=zoom)
    finally:
        os.remove(outfile)


@cmd.extend
def load_aln(filename, object=None, mobile=None, target=None, mobile_id=None,
             target_id=None, format='', transform=0, quiet=1, *, _self=cmd):
    '''
    DESCRIPTION
        Load a pairwise alignment from file and apply it to two loaded structures.
    USAGE
        load_aln filename [, object [, mobile [, target [, mobile_id [, target_id [, format ]]]]]]
    ARGUMENTS
        filename = string: alignment file
        object = string: name of the object {default: filename prefix}
        mobile, target = string: atom selections {default: ids from alignment file}
        mobile_id, target_id = string: ids from alignment file {default: first two}
        format = string: file format, see http://biopython.org/wiki/AlignIO
        {default: guess from first line in file}
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    EXAMPLE
        >>> fetch 1bz4 1cpr, async=0
        >>> super 1bz4 and guide, 1cpr and guide, object=aln1, window=5
        >>> save /tmp/super.aln, aln1
        >>> delete aln1
        >>> load_aln /tmp/super.aln, aln2, format=clustal
    '''
    import os
    from . import one_letter
    from .seqalign import needle_alignment, alignment_mapping, alignment_read

    quiet = int(quiet)
    if object is None:
        object = os.path.basename(filename).rsplit('.', 1)[0]

    # load alignment file
    alignment = alignment_read(cmd.exp_path(filename), format)
    aln_dict = dict((record.id, record) for record in alignment)
    mobile_record = alignment[0] if mobile_id is None else aln_dict[mobile_id]
    target_record = alignment[1] if target_id is None else aln_dict[target_id]

    # guess selections from sequence identifiers (if not given)
    if mobile is None:
        mobile = mobile_record.id
    if target is None:
        target = target_record.id

    try:
        mobile_obj = _self.get_object_list('(' + mobile + ')')[0]
        target_obj = _self.get_object_list('(' + target + ')')[0]
    except (
            IndexError,  # empty list (no atoms in selection)
            TypeError,  # None (Selector-Error)
            CmdException,  # future-proofing (Selector-Error)
    ):
        raise CmdException(
            f'selection "{mobile}" or "{target}" does not exist') from None

    # get structure models and sequences
    mobile_model = _self.get_model('(%s) and guide' % mobile)
    target_model = _self.get_model('(%s) and guide' % target)
    mobile_sequence = ''.join(one_letter.get(a.resn, 'X')
                              for a in mobile_model.atom)
    target_sequence = ''.join(one_letter.get(a.resn, 'X')
                              for a in target_model.atom)

    # align sequences from file to structures
    mobile_aln = needle_alignment(str(mobile_record.seq), mobile_sequence)
    target_aln = needle_alignment(str(target_record.seq), target_sequence)

    # get index mappings
    mobile_aln_i2j = dict(alignment_mapping(*mobile_aln))
    target_aln_i2j = dict(alignment_mapping(*target_aln))
    record_i2j = alignment_mapping(mobile_record, target_record)

    # build alignment list
    r = []
    for i, j, in record_i2j:
        if i in mobile_aln_i2j and j in target_aln_i2j:
            i = mobile_aln_i2j[i]
            j = target_aln_i2j[j]
            r.append([
                (mobile_obj, mobile_model.atom[i].index),
                (target_obj, target_model.atom[j].index),
            ])

    cmd.set_raw_alignment(object, r, int(transform))
    return r


# Register extensions
try:
    from pymol.importing import loadfunctions
    loadfunctions.setdefault("csv", load_csv)
    loadfunctions.setdefault("top", load_topol)
    loadfunctions.setdefault("ndx", load_ndx)
    loadfunctions.setdefault("aln", load_aln)
except ImportError:
    pass
