#!/usr/bin/env python3
from pymol import cmd, cgo
import numpy as np
from scipy import spatial
import collections
import re

def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1, color='blue red', name=''):      
  """
DESCRIPTION
  Create a CGO arrow between two picked atoms.
ARGUMENTS
  atom1 = string: single atom selection or list of 3 floats {default: pk1}
  atom2 = string: single atom selection or list of 3 floats {default: pk2}
  radius = float: arrow radius {default: 0.5}
  gap = float: gap between arrow tips and the two atoms {default: 0.0}
  hlength = float: length of head
  hradius = float: radius of head
  color = string: one or two color names {default: blue red}
  name = string: name of CGO object
  """
  from chempy import cpv

  radius, gap = float(radius), float(gap)
  hlength, hradius = float(hlength), float(hradius)

  try:
    color1, color2 = color.split()
  except:
    color1 = color2 = color
  color1 = list(cmd.get_color_tuple(color1))
  color2 = list(cmd.get_color_tuple(color2))

  def get_coord(v):
    if not isinstance(v, str):
      return v
    if v.startswith('['):
      return cmd.safe_list_eval(v)
    return cmd.get_atom_coords(v)

  xyz1 = get_coord(atom1)
  xyz2 = get_coord(atom2)
  normal = cpv.normalize(cpv.sub(xyz1, xyz2))

  if hlength < 0:
    hlength = radius * 3.0
  if hradius < 0:
    hradius = hlength * 0.6

  if gap:
    diff = cpv.scale(normal, gap)
    xyz1 = cpv.sub(xyz1, diff)
    xyz2 = cpv.add(xyz2, diff)

  xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

  obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
        [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
        [1.0, 0.0]

  if not name:
      name = cmd.get_unused_name('arrow')

  cmd.load_cgo(obj, name)

  return


def sum_partial_charges(selection="(all)", quiet=1, _self=cmd):
  _util_sum_pc = [0.0]
  _self.iterate(selection, "_util_sum_pc[0] += partial_charge", space=locals())
  result = _util_sum_pc[0]
  if not int(quiet):
    print(" util.sum_partial_charges: sum = %0.4f" % result)
  return result


def load_topol(file, selection="(all)", _self=cmd) -> None:
  """
DESCRIPTION
  Read partial_charges from GROMACS topology file.
USAGE
  load_topol file [, selection]
ARGUMENTS
  file = str: Input file name
  selection = str: Selection to load data to. (default: all)
PYTHON API
  cmd.load_topol(...) -> None
  """
  # Properties we want to import
  atom_list = collections.defaultdict(float)
  # Process topology file. We are mostly interested in `atoms` directive

  current_keyword = ""
  for line in open(file, "r"):
    # Remove comments after ';' symbol 
    line = line.split(";")[0].strip()
    if not line:
      continue
    def get_keyword(string):
      """ Scan a string for directive """
      k = re.search(r'\[(.*?)\]', string)
      return k.string.split()[1] if k else None
    # Groups are in square parentheses
    if "[" in line:
      current_keyword = get_keyword(line)

    elif current_keyword == "atoms":
      # nr, type, resnr, residue, atom, cgnr, pcharge, mass
      col = line.split()
      try:
        index, pcharge = int(col[0]), float(col[6])
        atom_list[index] = pcharge
      except:
        print(" Util: Cannot parse line: '%s'" % line)

  _self.alter(selection, "partial_charge = atom_list[index]", space=locals())
  sum_partial_charges(selection, quiet=0, _self=_self)

  return


def get_longest(selection="(all)", vis=False, _self=cmd) -> float:
  """ Get longest distance in molecule """
  xyz = np.array(_self.get_coords(selection, 1))

  # Find 2 points with largest distance - Convex Hull
  xyz_hull = xyz[spatial.ConvexHull(xyz).vertices]
  dist_mat = spatial.distance_matrix(xyz_hull, xyz_hull)
  i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
  x = xyz_hull[i] - xyz_hull[j]

  if vis:
    cgo_arrow(xyz_hull[i].tolist(), xyz_hull[j].tolist(), color="green")

  # Old fashined way
  # num = len(xyz)
  # r_max = 0.0
  # for i in range(0, num-1):
  #   for j in range(i+1, num):
  #     r = np.linalg.norm(xyz[i] - xyz[j])
  #     if r > r_max:
  #       r_max = r
  #       ni, nj = i, j
  # x = xyz[ni] - xyz[nj]

  return x 


def get_dipole(selection="(all)", vis=False, _self=cmd) -> float:
  """ Calculate dipole moment vector. """
  com = np.array(_self.centerofmass(selection))
  xyz = np.array(_self.get_coords(selection, 1))

  pcharge = []
  _self.iterate(selection, "pcharge.append(partial_charge)", space=locals())

  dipole = np.array([0., 0., 0.])
  for i in range(len(xyz)):
    dipole += pcharge[i] * (xyz[i] - com)

  # Transform units [e*nm] to [D]
  dipole /= 0.2081943 

  if vis:
    cgo_arrow(com.tolist(), dipole.tolist())

  return dipole


def get_angle(u, v) -> float:
  """ Get angle in degrees between vector u and v """
  r = np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)
  return np.degrees(np.arccos(r))


# ----------------------------------------------

# 1JKD w/ oplsaa force-field has dipole moment:
#  Total < M_x > = 150.078 Debye
#  Total < M_y > = 30.0526 Debye
#  Total < M_z > = -79.0612 Debye

sele = "conf"

# Clean up
cmd.delete("arrow01")
cmd.delete("arrow02")

# Load structure and topology
cmd.load("conf.gro")
load_topol("topol.top")

# Get longest distance vector
x = get_longest(sele)
print(x)

d = get_dipole(sele)
print(d)

result = get_angle(x, d)
result = min(result, abs(result - 180))
print(result)
