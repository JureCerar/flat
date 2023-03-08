#!/usr/bin/env python3
from pymol import cmd, CmdException

def load_or_update(molstr, name, sele, state, _self):
  with _self.lockcm:
    _load_or_update(molstr, name, sele, state, _self)


def _load_or_update(molstr, name, sele, state, _self):
  update = not name
  if update:
    name = _self.get_unused_name('_minimized')
  else:
    _self.delete(name)
  _self.load_raw(molstr, 'mol', name, 1, zoom=0)
  try:
    from psico.fitting import xfit
    xfit(name, sele, 1, state, match='none', cycles=100, guide=0)
  except ImportError:
    _self.fit(name, sele, 1, state, cycles=5, matchmaker=-1)
  except Exception as e:
    print('xfit failed: {}'.format(e))
  if update:
    _self.update(sele, name, state, 1, matchmaker=0)
    _self.delete(name)
  return


def get_fixed_indices(selection, state, _self):
  fixed_list = []
  _self.iterate_state(
    state, selection,
    '_append(flags & 0x8)',
    space={'_append': fixed_list.append}
  )
  return [idx for (idx, fixed) in enumerate(fixed_list) if fixed]


def randomize_coords_if_collapsed(selection, state, fancy=True, _self=cmd):
  '''
  If all coordinates are the same (collapsed into one point), then
  randomize them.
  :param fancy: Arrange atoms in a circle (this works better for openbabel)
  :type fancy: bool
  '''
  import numpy.random
  coords = _self.get_coords(selection, state)
  if len(coords) < 2 or coords.std(0).sum() > 1e-3:
    return
  if fancy:
    # puts x,y on a circle
    angles = numpy.linspace(0, 2 * numpy.pi, len(coords), False)
    width = len(coords)**(1 / 3.)
    coords[:, 0] += numpy.sin(angles) * width
    coords[:, 1] += numpy.cos(angles) * width
  coords += numpy.random.random_sample(coords.shape) - 0.5
  _self.load_coords(coords, selection, state)
  return


@cmd.extend
def minimize_ob(selection='enabled', state=-1, ff='UFF', nsteps=500,
                conv=0.0001, cutoff=0, cut_vdw=6.0, cut_elec=8.0,
                name='', quiet=1, _self=cmd):
  '''
DESCRIPTION
  Emergy minimization with openbabel
  Supports fixed atoms (flag fix)
ARGUMENTS
  selection = str: atom selection
  state = int: object state {default: -1}
  ff = GAFF|MMFF94s|MMFF94|UFF|Ghemical: force field {default: UFF}
  nsteps = int: number of steps {default: 500}
  '''
  import openbabel as ob
  try:
    # OB 3.x
    from openbabel import openbabel as ob
  except ImportError:
    # OB 2.x
    pass

  state = int(state)
  sele = _self.get_unused_name('_sele')
  natoms = _self.select(sele, selection, 0)
  try:
    if natoms == 0:
      raise CmdException('empty selection')
    randomize_coords_if_collapsed(sele, state, _self=_self)
    ioformat = 'mol'
    molstr = _self.get_str(ioformat, sele, state)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats(ioformat, ioformat)
    mol = ob.OBMol()
    obconversion.ReadString(mol, molstr)
    # add hydrogens
    orig_ids = [a.GetId() for a in ob.OBMolAtomIter(mol)]
    mol.AddHydrogens()
    added_ids = set(
      a.GetId() for a in ob.OBMolAtomIter(mol)
    ).difference(orig_ids)
    consttrains = ob.OBFFConstraints()
    consttrains.Setup(mol)
    # atoms with "flag fix"
    fixed_indices = get_fixed_indices(sele, state, _self)
    for idx in fixed_indices:
      consttrains.AddAtomConstraint(idx + 1)
    # setup forcefield (one of: GAFF, MMFF94s, MMFF94, UFF, Ghemical)
    ff = ob.OBForceField.FindForceField(ff)
    if ff is None:
      raise CmdException(
        "FindForceField returned None, please check BABEL_LIBDIR and BABEL_DATADIR"
      )
    ff.Setup(mol, consttrains)
    if int(cutoff):
      ff.EnableCutOff(True)
      ff.SetVDWCutOff(float(cut_vdw))
      ff.SetElectrostaticCutOff(float(cut_elec))
    # run minimization
    ff.SteepestDescent(int(nsteps) // 2, float(conv))
    ff.ConjugateGradients(int(nsteps) // 2, float(conv))
    ff.GetCoordinates(mol)
    # remove previously added hydrogens
    for hydro_id in added_ids:
      mol.DeleteAtom(mol.GetAtomById(hydro_id))
    molstr = obconversion.WriteString(mol)
    load_or_update(molstr, name, sele, state, _self)
    if not int(quiet):
      print(' Energy: %8.2f %s' % (ff.Energy(), ff.GetUnit()))
  finally:
    _self.delete(sele)
  return


# minimize_ob(selection='sele', state=-1, ff='UFF', nsteps=500,
#             conv=0.0001, cutoff=0, cut_vdw=6.0, cut_elec=8.0,
#             name='', quiet=1, _self=cmd)
