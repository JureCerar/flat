from pymol import cmd, CmdException

@cmd.extend
def renumber(selection, start=1, quiet=False):
  """
DESCRIPTION
  Set residue numbering (resi) based on connectivity.
ARGUMENTS
  selection = string: atom selection to renumber {default: all}
  start = integer: counting start {default: 1}
  """
  import sys
  start, quiet = int(start), bool(quiet)
  model = cmd.get_model(selection)
  limit = sys.getrecursionlimit()
  sys.setrecursionlimit(10000)
  cmd.iterate(
    selection,
    'next(atom_it).model = model',
    space={
      'atom_it': iter(model.atom),
      'next': next
    }
  )
  startatom = model.atom[0]
  for atom in model.atom:
    atom.adjacent = []
    atom.visited = False
  for bond in model.bond:
    atoms = [model.atom[i] for i in bond.index]
    atoms[0].adjacent.append(atoms[1])
    atoms[1].adjacent.append(atoms[0])
  minmax = [start, start]

  def traverse(atom, resi):
    atom.resi = resi
    atom.visited = True
    for other in atom.adjacent:
      if other.visited:
        continue
      if (atom.name, other.name) in [('C', 'N'), ("O3'", 'P')]:
        minmax[1] = resi + 1
        traverse(other, resi + 1)
      elif (atom.name, other.name) in [('N', 'C'), ('P', "O3'")]:
        minmax[0] = resi - 1
        traverse(other, resi - 1)
      elif (atom.name, other.name) not in [('SG', 'SG')]:
        traverse(other, resi)
    return
  
  traverse(startatom, start)
  cmd.alter(
    selection,
    'resi = next(atom_it).resi',
    space={
      'atom_it': iter(model.atom),
      'next': next
    }
  )
  sys.setrecursionlimit(limit)
  if not quiet:
    print(' Renumber: range (%d to %d)' % tuple(minmax))
  return tuple(minmax)


@cmd.extend
def graft (mobile, target, output=None):
  """
DESCRIPTION
  Graft mobile structure on target structure.
ARGUMENTS
  mobile = string: atom selection of mobile section 
  target = string: atom selection of target section 
  output = string: output structure name {default: None}
  """
  if len(cmd.get_chains(mobile)) != 1:
    raise CmdException("Mobile selection must contain only one chain")
  if len(cmd.get_chains(target)) != 1:
    raise CmdException("Target selection must contain only one chain")

  _mobile = cmd.get_unused_name("mobile")
  cmd.copy_to(_mobile, mobile, rename="", zoom=0)
  mobile_len = len(cmd.get_model(_mobile).get_residues())

  _target = cmd.get_unused_name("target")
  cmd.copy_to(_target, target, rename="", zoom=0)
  target_len = len(cmd.get_model(_target).get_residues())

  # Clean-up mobile section to match target
  model = cmd.get_model(_target)
  chain, segi = model.atom[1].chain, model.atom[1].segi
  cmd.alter(_mobile, f"chain='{chain}'; segi='{segi}'")
  cmd.remove(f"{_mobile} & name OXT")

  # cealign metod works best for aligment
  cmd.extra_fit(_mobile, _target, "cealign", mobile_state=-1, target_state=-1)
  align = cmd.get_unused_name("align")
  score = cmd.align(_mobile, _target, quiet=0, object=align)
  if score[0] > 100.0:
    raise Warning("RMSD value is high. Are you sure this is correct aligment?")

  # Renumber target and mobile to correct numbering
  model = cmd.get_model(f"{align} & {_target}")
  first, last = int(model.atom[0].resi), int(model.atom[-1].resi)
  renumber(_target, start=1)
  renumber(_mobile, start=first)
  cmd.remove(f"{_target} & resi {first}-{last}")
  cmd.delete(align)

  # Combine to new objects
  if not output:
    output = cmd.get_unused_name("out")
  cmd.copy_to(output, _target, rename="", zoom=0)
  cmd.copy_to(output, _mobile, rename="", zoom=0)
  cmd.delete(_target)
  cmd.delete(_mobile)

  # Bond mobile and target chains.
  try:
    cmd.bond(f"/{output}///{last}/C", f"/{output}///{last+1}/N")
    cmd.bond(f"/{output}///{first-1}/C", f"/{output}///{first}/N")
  except:
    pass

  # Minimize around new bond to avoid any clashes.
  cmd.rebuild(output)
  # try:
  #   from minimize import minimize_ob
  #   sele = cmd.get_unused_name("sele")
  #   cmd.select(sele, f"{output} & resi {last} around {10.0}")
  #   minimize_ob(sele, state=-1, ff="UFF", nsteps=500, quiet=0)
  # finally:
  #   cmd.delete(sele)

  return


cmd.delete("all")
mobile = "1gig"
target = "1hzh"
cmd.load("conf/{}.pdb".format(mobile), mobile)
cmd.load("conf/{}.pdb".format(target), target)


graft(f"{mobile} & chain H", f"{target} & chain H")
graft(f"{mobile} & chain L", f"{target} & chain L")

cmd.delete(mobile)
cmd.delete(target)

cmd.save("out.pdb")

# cmd.delete(mobile)
# cmd.delete(target)
