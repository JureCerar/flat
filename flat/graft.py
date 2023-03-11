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
      'next': next,
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
      'next': next,
    }
  )
  sys.setrecursionlimit(limit)
  if not quiet:
    print(' Renumber: range (%d to %d)' % tuple(minmax))
  return tuple(minmax)


@cmd.extend
def graft (mobile, target, out=None, minimize=True):
  """
DESCRIPTION
  Graft mobile structure on target structure. For successfulk graft,
  similarity between aligned mobile and target structure must be > 80 %.
ARGUMENTS
  mobile = string: atom selection of mobile section 
  target = string: atom selection of target section 
  out = string: output structure name {default: None}
  minimize = bool: Do a short minimization of grafted structure {default: True}
  """
  minimize = bool(minimize)
  if minimize:
    from minimize import minimize_ob

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

  # Clean-up mobile section; chain and segi must match target section.
  model = cmd.get_model(_target)
  chain, segi = model.atom[0].chain, model.atom[0].segi
  cmd.alter(_mobile, f"chain='{chain}'; segi='{segi}'")
  cmd.remove(f"{_mobile} & name OXT")

  # 'cealign' method works best for aligment
  cmd.extra_fit(_mobile, _target, "cealign", mobile_state=-1, target_state=-1)
  align = cmd.get_unused_name("align")
  score = cmd.align(_mobile, _target, quiet=0, object=align)
  if score[0] > 100.0:
    raise Warning("RMSD value is high. Are you sure this is correct aligment?")

  # Renumber target and mobile to correct numbering
  model = cmd.get_model(f"{align} & {_target}")
  first, last = int(model.atom[0].resi), int(model.atom[-1].resi)

  if not out:
    out = cmd.get_unused_name("out")
  renumber(_mobile, start=first)
  cmd.copy_to(out, _mobile, rename="", zoom=0)

  # Split target section into left and right objects.
  # Correctly renumber objects and then combine/join them.
  # Do a energy minimization around newly bonded sections.
  if first > 1:
    target_left = cmd.get_unused_name("target_left")
    cmd.select(target_left, f"{_target} & resi 1-{first-1}")
    renumber(target_left, start=1)
    cmd.copy_to(out, target_left, rename="", zoom=0)
    cmd.bond(
      f"/{out}///{first-1}/C",
      f"/{out}///{first}/N",
    )
    cmd.delete(target_left)
    if minimize:
      minimize_ob(
        f"{out} & resi {first} around {10.0}",
        ff="UFF",
        conv=0.00001,
        nsteps=500,
        quiet=0,
      )

  if last < target_len:
    target_right = cmd.get_unused_name("target_right")
    cmd.select(target_right, f"{_target} & resi {last+1}-{target_len}")
    renumber(target_right, start=mobile_len+first)
    cmd.copy_to(out, target_right, rename="", zoom=0)
    cmd.bond(
      f"/{out}///{first+mobile_len-1}/C",
      f"/{out}///{first+mobile_len}/N",
    )
    cmd.delete(target_right)
    if minimize:
      minimize_ob(
        f"{out} & resi {first+mobile_len} around {10.0}",
        ff="UFF",
        conv=0.00001,
        nsteps=500,
        quiet=0,
      )

  # Clean-up objects (in order they were created)
  cmd.delete(_mobile)
  cmd.delete(_target)
  cmd.delete(align)

  return


cmd.delete("all")
mobile = "1dlf"
target = "1gig"

cmd.load("conf/{}.pdb".format(mobile), mobile)
cmd.load("conf/{}.pdb".format(target), target)

graft(f"{mobile} & chain H", f"{target} & chain H", minimize=False)
graft(f"{mobile} & chain L", f"{target} & chain L", minimize=False)

cmd.delete(mobile)
cmd.delete(target)

cmd.save("out.pdb")

