#!/usr/bin/env python3
from pymol import cmd
import collections

@cmd.extend
def saveBfact(file=None, selection="all", var="b", _self=cmd) -> None:
  """
DESCRIPTION
USAGE
ARGUMENTS
NOTES
PYTHON API
SEE ALSO
  """
  bfact = collections.defaultdict(float)
  _self.iterate(selection, "bfact[model,segi,chain,resn,resi,name] =" + var, space=locals())

  with open(file, "w") as handle:
    print("model", "segi", "chain", "resn", "resi", "name", "bfact", sep=",", file=handle)
    for key, val in bfact.items():
      print(",".join(key), val, sep=",", file=handle)

  return


@cmd.extend
def loadBfact(file, selection="all", var="b", vis=1, _self=cmd):
  """
DESCRIPTION
USAGE
ARGUMENTS
NOTES
PYTHON API
SEE ALSO
  """
  vis = bool(vis)
  bfact = collections.defaultdict(float)
  with open(file, "r") as handle:
    line = handle.readline()
    for line in handle:
      col = line.split(",")
      key, val = tuple(col[0:6]), float(col[6])
      bfact[key] = val

  _self.alter(selection, var + "= bfact[model,segi,chain,resn,resi,name]", space=locals())

  if vis:
    obj = _self.get_object_list(selection)[0]
    ramp = _self.get_unused_name("ramp")
    _self.spectrum(var, "rainbow", selection)
    _self.ramp_new(ramp, obj, range=[min(bfact.values()), max(bfact.values())], color="rainbow")

  return