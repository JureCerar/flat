from pymol import cmd

@cmd.extend
def chi (selection="all", state=1, quiet=False) -> dict:
  """
DESCRIPTION
  Get chi side-chain dihedral angle for each residue.
ARGUMENTS
  selection = str: Atom selection {default: all}
  state = int: Object state (0 for all states). {default: 0}
  quiet = bool: Print results to log window. {default: True}
PYTHON API
  cmd.chi(...) -> dict
  """
  chi_angle = {}
  atoms = ["N", "CA", "CB", "CG"]
  for obj in cmd.get_names('public_nongroup_objects', selection=selection):
    for res in cmd.get_model(f"bca. {obj}").atom:
      key = (obj, res.segi, res.chain, res.resi)
      residue = "/{}/{}/{}/{}/".format(*key)
      try:
        angle = cmd.get_dihedral(
          residue + atoms[0],
          residue + atoms[1],
          residue + atoms[2],
          residue + atoms[3],
          state=state
        )
        chi_angle[key] = angle
      except:
        continue
  if not quiet:
    for key, val in chi_angle.items():
      print("/{}/{}/{}/{}/ {:.3f}".format(*key, val))
  return chi_angle