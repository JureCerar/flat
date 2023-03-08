#!/usr/bin/env python3
from pymol import cmd, util

@cmd.extend
def beautify(selection="all") -> None:
  """
  DESCRIPTION
    Display molecules in pretty representation (as name implies).
  ARGUMENTS
    selection = str: Atom selection. {default: all}
  """
  cmd.hide("everything")
  # Show protein
  sele = cmd.get_unused_name("_sele")
  cmd.select(sele, "%s and polymer.protein()" % selection)
  cmd.show("cartoon", sele)
  cmd.show("licorice", "%s & (sidechain | name CA) & ! hydrogen" % sele)
  # Show not-protein
  cmd.show("wire", "! %s" % sele)
  # Color atoms
  util.cbc()
  cmd.color("red", "element O")
  cmd.color("blue", "element N")
  cmd.color("yellow", "element S")
  cmd.delete(sele)
  return None