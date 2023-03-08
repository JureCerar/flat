#!/usr/bin/env python3

def read_ndx(file) -> dict:
  """
  DESCRIPTION
    Read a GROMACS index file as dictionary.
  USAGE
    read_ndx file
  ARGUMENTS
    file = str: Input file name.
  PYTHON API
    cmd.read_ndx(...) -> dict
  """
  import re
  group_name = None
  groups = dict()
  for line in open(file, "r"):
    # Remove comments if any
    line = line.split(";")[0]
    # Scan a string for directive
    def get_group(string):
      k = re.search(r'\[(.*?)\]', string)
      return k.string.split()[1] if k else None
    # Process lines
    if "[" in line:
      group_name = get_group(line)
      groups[group_name] = []
    elif group_name:
      groups[group_name] += [int(i) for i in line.split()]
  return groups

def write_ndx(groups, file=None) -> None:
  """
  DESCRIPTION
    Write dictionary in into a GROMACS index file.
  USAGE
    write_ndx groups [, file ]
  ARGUMENTS
    groups = dict: Dictionary containing atom indeces
    file = str: Output file name. (default: STDOUT)
  PYTHON API
    cmd.write_ndx(...) -> None
  """
  handle = open(file, "w") if file else None

  nlines = 15 # Number of lines
  for name, group in groups.items():
    print("[ %s ]" % name, file=handle)
    for i in range(0, len(group), nlines):
      print(" ".join("{:5}".format(x) for x in group[i:i+nlines]), file=handle)

  if file:
    handle.close()

  return

if __name__ == "pymol" :
  from pymol import cmd
    
  def load_ndx(file, _self=cmd) -> None:
    """
    DESCRIPTION
      Create selections from a GROMACS index file.
    USAGE
      load_ndx file
    ARGUMENTS
      file = str: Input file name
    PYTHON API
      cmd.load_ndx(...) -> None
    SEE ALSO
      save_ndx
    """
    groups = read_ndx(file)

    for name, group in groups.items():
      # Pymol does not like '&' and '|' characters and reserved keywords
      name = name.replace("&", "_and_")
      name = name.replace("|", "_or_")
      if name.lower() in ["sidechain", "backbone"]:
        name = name + "_"

      if len(group) > 0:
        print(" Loading group: '%s' (%s atoms)" % (name, len(group)))
        # Pymol does not like long selection strings so we select the group iteratively
        every = 15
        _self.select("(%s)" % name, "id %i" % group[0])
        for i in range(1, len(group), every):
          buffer = "index " + "+".join(str(x) for x in group[i:i+every] )
          _self.select(name, "(%s) | %s" % (name, buffer), quiet=1)

      else:
        print(" Empty group: '%s'" % name)

    return

  def save_ndx(file, _self=cmd) -> None:
    """
    DESCRIPTION
      Save all the selections into a GROMACS index file.
    USAGE
      save_ndx file
    ARGUMENTS
      file = str: Output file name
    PYTHON API
      cmd.save_ndx(...) -> None
    SEE ALSO
      load_ndx
    """
    names = _self.get_names("public_selections")
    groups = {}

    for name in names:
      indices = []
      _self.iterate(name, "indeces.append(index)", space=locals())
      if indices:
        groups[name] = indices

    write_ndx(groups, file)

    return

  # Declare the new commands in Pymol
  cmd.extend("load_ndx", load_ndx)
  cmd.extend("save_ndx", save_ndx)
