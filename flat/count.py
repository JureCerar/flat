#!/usr/bin/env python3
from pymol import cmd

@cmd.extend
def count(selection="all") -> None:
  """
  DESCRIPTION
    Count number of atoms, residues, and mass in selection.
  USAGE
    count [ selection ]
  ARGUMENTS
    selection = str: Atom selection. {default: all}
  """
  # Get selection
  sele = cmd.get_unused_name("_sele")
  cmd.select(sele, selection)
  # Get properties
  model = cmd.get_model(selection)
  print("Residue(s): %s" % len(model.get_residues()))
  print("Atom(s): %s" % model.nAtom)
  mass = model.get_mass()
  if mass < 1000:
    print("Mass: %.3f Da" % mass)
  else:
    print("Mass: %.3f kDa" % (mass/1000))
  cmd.delete(sele)
  return None
  
  

  