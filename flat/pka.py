from pymol import cmd, CmdException
import propka.run as pk
import tempfile
import os

# Ref: https://propka.readthedocs.io/en/latest/api.html
# Ref: https://tuilab.github.io/tui/update/2017/02/21/propka.html

# To-Do:
# - [ ] Iterate over states & average values
# - [ ] Return dict() or create Residue() class
# - [ ] Plot titration curve and buffer capacity
# - [ ] Get isoelectric poin
# 

def get_propka(selection):
  """ Calculate pKa """

  filename = os.path.join(tempfile.mkdtemp(), f'something.pdb')
  optargs = []
  
  # print(filename)
  cmd.save(filename, selection, state=1)

  protein = pk.single(filename, optargs, None, False)

  cmd.label()

  residues_out = ["ASP", "GLU", "HIS", "TYR", "LYS", "ARG", "N+", "C-"]

  # Loop over atoms
  conformation = "AVR"
  for group in protein.conformations[conformation].groups:
    if group.residue_type in residues_out:

      resi = group.atom.res_num
      resn = group.atom.res_name
      chain = group.atom.chain_id
      name = group.atom.name
      buried = group.pka_value
      pka_value = group.pka_value
      pka_model = group.model_pka

      label = f"'pKa={pka_value:.2f}({pka_model-pka_value:.2f}) Bu:{buried:.0f}%'"
      selection = f"///{chain}/{resn}`{resi}/{name}"

      # print(selection)
      cmd.label(selection, label)
      cmd.show("spheres", selection)

  cmd.set("label_size", 6)
  cmd.set("sphere_scale", 0.25, "(all)")


# get_propka("1jkd")
# print(group.__dict__)