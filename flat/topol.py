from pymol import cmd, util
import collections
import re

@cmd.extend
def load_topol(file, selection="(all)", _self=cmd):
  """ Comment """
  # Properties we want to import
  atom_list = collections.defaultdict(float)
  # Process topology file. We are mostly interested in `atoms` directive
  current_keyword = ""
  for line in open(file, "r"):
    line = line.split(";")[0].strip()
    if not line:
      continue
    # Scan a string for directive
    def get_keyword(string):
      k = re.search(r'\[(.*?)\]', string)
      return k.string.split()[1] if k else None
    # Update current section keyword
    keyword = get_keyword(line)
    if keyword:
      current_keyword = keyword
      continue
    if current_keyword == "atoms":
      # nr, type, resnr, residue, atom, cgnr, pcharge, mass
      col = line.split()
      try:
        index, pcharge = int(col[0]), float(col[6])
        atom_list[index] = pcharge
      except:
        print("Cannot parse line: '%s'" % line)
  # Set atom properties
  _self.alter(selection, "partial_charge = atom_list[index]", space=locals())
  util.sum_partial_charges(selection, quiet=0, _self=_self)
  return