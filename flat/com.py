#!/usr/bin/env python3
from pymol import cmd

@cmd.extend
def get_com(selection, state=1, quiet=0) -> float:
  """
  Get centre-of-mass (COM) for selection and state.
  """
  state = int(state)
  quiet = bool(quiet)
  
  model = cmd.get_model(selection, state)
  total_mass = 0.
  x = 0.
  y = 0.
  z = 0. 
  
  for atom in model.atom:
    mass = atom.get_mass()
    x += atom.coord[0] * mass
    y += atom.coord[1] * mass
    z += atom.coord[2] * mass
    total_mass += mass
 
  if total_mass: 
    x = x / total_mass
    y = y / total_mass
    z = z / total_mass
  
  if not quiet:
    print("COM[%i]={%f,%f,%f}" % (state, x, y, z))
  
  return x, y, z


@cmd.extend
def show_com(selection="all", state=None, object=None, **kwargs) -> None:
  """
  DESCRIPTION
    Display centre-of-mass (COM) for selected atoms and state.
  ARGUMENTS
    selection = str: Atom selection. {default: all}
    state = int: Selection state. {default: None}
    object = str: Output object name. {default: None}
    kwargs = *: Optional parameters passed to 'cmd.pesudoatom()'.
  """
  state = int(state)
  
  # Get selection
  sele = cmd.get_unused_name("_sele")
  cmd.select(sele, selection)
  
  # Generate object name
  if object == None:
    try:
      object = cmd.get_legal_name(sele)
      object = cmd.get_unused_name(object + "_COM", 0)
    except:
      object = "COM"
  cmd.delete(object)

  # Get COM for defined state  
  if state != None:
    x, y, z = get_com(sele, quiet=1)
    cmd.pseudoatom(object, pos=[x, y, z], vdw=0.25, **kwargs)
    cmd.show("spheres", object)
  else:
    for state in range(cmd.count_states()):
      x, y, z = get_com(sele, state=state+1, quiet=1)
      cmd.pseudoatom(object, pos=[x, y, z], vdw=0.25, **kwargs)
      cmd.show("spheres", object)
  
  return None