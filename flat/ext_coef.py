#!/usr/bin/env python3
from pymol import cmd, util

@cmd.extend
def ext_coef(selection="all") -> None:
  """
  DESCRIPTION
    Calculate extinction coefficient and absorbance for native and folded protein.
  ARGUMENTS
    selection = str: Atom selection. {default: all}
  LITERATURE
    [1] C.N. Pace, et. al, Protein Sci., 1995, DOI: https://doi.org/10.1002/pro.5560041120 
    [2] H. Edelhoch, Biochemistry, 1967, DOI: https://doi.org/10.1021/bi00859a010
    [3] C.G. Gill & P.H. von Hippel, Anal. Biochem., 1989, DOI: https://doi.org/10.1016/0003-2697(89)90602-7
  """
  # count number of residues in model
  tyr = len(cmd.get_model("({}) & resn TYR".format(selection)).get_residues())
  trp = len(cmd.get_model("({}) & resn TRP".format(selection)).get_residues())
  cys = len(cmd.get_model("({}) & resn CYS+CYS2".format(selection)).get_residues())
  cys2 = cmd.get_model("({}) & (CYS+CYS2/SG and bound_to CYS+CYS2/SG)".format(selection)).nAtom
  
  # If not higher structure data is avaliable (i.e. calculating from fasta)
  cys2 = cys if cys2 == 0 else cys2
  
  # Get protein mass
  mass = cmd.get_model(selection).get_mass()
  
  # Calculate extinction coefficeint and absorbance at 1 mg/ml for native protein
  coef = trp * 5500 + tyr * 1490 + cys2 * 62.5
  abs = coef / mass
  
  # Calculate for denatured protein
  coef0 = trp * 5500 + tyr * 1490
  abs0 = coef0 / mass
  
  # Return results in nice format
  print(
    "Extinction coefficients are in units of M^-1 cm^-1, at 280 nm measured in water.",
    "",
    "Ext. coefficeint: {:.0f}".format(coef),
    "Abs 0.1% (=1 mg/mL): {:.3f}, for native protein.".format(abs),
    "",
    "Ext. coefficeint: {:.0f}".format(coef0),
    "Abs 0.1% (=1 mg/mL): {:.3f}, for denatured protein.".format(abs0),
    sep="\n"
  )
    
  return None
