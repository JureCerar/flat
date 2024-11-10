# Copyright (C) 2023-2024 Jure Cerar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from pymol import cmd, CmdException
import io


def get_fixed_indices(selection, state, _self):
    """Fix pymol idices"""
    fixed_list = []
    _self.iterate_state(
        state,
        selection,
        "_append(flags & 0x8)",
        space={"_append": fixed_list.append}
    )
    return [idx for (idx, fixed) in enumerate(fixed_list) if fixed]


def load_or_update(molstr, name, fmt, sele, state, _self):
    """Load or update structure of existing object"""
    with _self.lockcm:
        update = not name

        if update:
            name = _self.get_unused_name("_minimized")
        else:
            _self.delete(name)

        _self.load_raw(molstr, fmt, name, 1, zoom=0)
        _self.fit(name, sele, 1, state, cycles=5, matchmaker=-1)

        if update:
            _self.update(sele, name, state, 1, matchmaker=0)
            _self.delete(name)


def randomize_coords_if_collapsed(selection, state, fancy=True, _self=cmd):
    """
    DESCRIPTION
        If all coordinates are the same (collapsed into one point), then
        randomize them.
    USAGE
        randomize_coords_if_collapsed [ selection [, state [, fancy ]]]
    ARGUMENTS
        selection = str: atom selection
        state = int: object state {default: -1}
        fancy = bool: Arrange atoms in a circle (this works better for openbabel)
    SOURCE
        From PSICO (c) Thomas Holder, Schrodinger Inc.
    """
    coords = _self.get_coords(selection, state)

    if len(coords) < 2 or coords.std(0).sum() > 1e-3:
        return

    import numpy.random

    if fancy:
        # puts x,y on a circle
        angles = numpy.linspace(0, 2 * numpy.pi, len(coords), False)
        width = len(coords)**(1 / 3.)
        coords[:, 0] += numpy.sin(angles) * width
        coords[:, 1] += numpy.cos(angles) * width

    coords += numpy.random.random_sample(coords.shape) - 0.5

    _self.load_coords(coords, selection, state)


@cmd.extend
def minimize_ob(selection="enabled", state=-1, ff="UFF", nsteps=500,
                conv=0.0001, cutoff=0, cut_vdw=6.0, cut_elec=8.0,
                name="", quiet=1, _self=cmd):
    """
    DESCRIPTION
        Energy minimization with OpenBabel. Supports fixed atoms (flag fix)
    USAGE
        minimize_ob [selection [, state [, ff [, nsteps [, ... ]]]]] 
    ARGUMENTS
        selection = str: atom selection
        state = int: object state {default: -1}
        ff = GAFF|MMFF94s|MMFF94|UFF|Ghemical: force field {default: UFF}
        nsteps = int: number of steps {default: 500}
    SOURCE
        From PSICO (c) Thomas Holder, Schrodinger Inc.
    """
    from openbabel import openbabel as ob  

    state = int(state)

    sele = _self.get_unused_name("_sele")
    natoms = _self.select(sele, selection, 0)

    try:
        if natoms == 0:
            raise CmdException("empty selection")

        randomize_coords_if_collapsed(sele, state, _self=_self)

        ioformat = "mol"
        molstr = _self.get_str(ioformat, sele, state)

        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats(ioformat, ioformat)

        mol = ob.OBMol()
        obconversion.ReadString(mol, molstr)

        # add hydrogens
        orig_ids = [a.GetId() for a in ob.OBMolAtomIter(mol)]
        mol.AddHydrogens()
        added_ids = set(a.GetId()
                        for a in ob.OBMolAtomIter(mol)).difference(orig_ids)

        consttrains = ob.OBFFConstraints()
        consttrains.Setup(mol)

        # atoms with "flag fix"
        fixed_indices = get_fixed_indices(sele, state, _self)
        for idx in fixed_indices:
            consttrains.AddAtomConstraint(idx + 1)

        # setup forcefield (one of: GAFF, MMFF94s, MMFF94, UFF, Ghemical)
        ff = ob.OBForceField.FindForceField(ff)
        if ff is None:
            raise CmdException("FindForceField returned None, please check "
                               "BABEL_LIBDIR and BABEL_DATADIR")
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
        load_or_update(molstr, name, "mol", sele, state, _self)

        if not int(quiet):
            print(" Energy: %8.2f %s" % (ff.Energy(), ff.GetUnit()))
    finally:
        _self.delete(sele)


@cmd.extend
def minimize_rdkit(selection="enabled", state=-1, ff="MMFF94", nsteps=500,
                   name="", quiet=1, _self=cmd):
    """
    DESCRIPTION
        Energy minimization with RDKit. Supports fixed atoms (flag fix)
    USAGE
        minimize_ob [selection [, state [, ff [, nsteps [, ... ]]]]] 
    ARGUMENTS
        selection = str: atom selection
        state = int: object state {default: -1}
        ff = MMFF94s|MMFF94|UFF: force field {default: MMFF94}
        nsteps = int: number of steps {default: 200}
    SOURCE
        From PSICO (c) Thomas Holder, Schrodinger Inc.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    state = int(state)

    sele = _self.get_unused_name("_sele")
    natoms = _self.select(sele, selection, 0)

    try:
        if natoms == 0:
            raise CmdException("empty selection")

        randomize_coords_if_collapsed(sele, state, _self=_self)

        molstr = _self.get_str("mol", sele, state)
        mol = Chem.MolFromMolBlock(molstr, True, False)

        if mol is None:
            raise CmdException("Failed to load molecule into RDKit. "
                               "Please check bond orders and formal charges.")

        # setup forcefield
        if ff.startswith("MMFF"):
            ff = AllChem.MMFFGetMoleculeForceField(
                mol,
                AllChem.MMFFGetMoleculeProperties(mol, ff, 0 if int(quiet) else 1)
            )
        elif ff == "UFF":
            ff = AllChem.UFFGetMoleculeForceField(mol)
        else:
            raise CmdException("unknown forcefield: " + ff)

        if ff is None:
            raise CmdException("forcefield setup failed")

        # atoms with "flag fix"
        for idx in get_fixed_indices(sele, state, _self):
            ff.AddFixedPoint(idx)

        # run minimization
        if ff.Minimize(int(nsteps)) != 0:
            print(" Warning: minimization did not converge")

        molstr = Chem.MolToMolBlock(mol)
        load_or_update(molstr, name, "mol", sele, state, _self)

        if not int(quiet):
            print(" Energy: %8.2f %s" % (ff.CalcEnergy(), "kcal/mol"))
    finally:
        _self.delete(sele)


@cmd.extend
def minimize_mm(selection="enabled", state=-1, ff="amber14-all.xml", nsteps=500,
                   name="", quiet=1, _self=cmd):
    """
    DESCRIPTION
        Energy minimization with OpenMM. Supports fixed atoms (flag fix)
    USAGE
        minimize [selection [, state [, ff [, nsteps [, ... ]]]]] 
    ARGUMENTS
        selection = str: atom selection {default: enabled}
        state = int: object state {default: -1}
        ff = str: OpenMM force field {default: amber14-all.xml}
        nsteps = int: number of steps {default: 200}
    """
    import openmm
    import openmm.app

    state, nsteps = int(state), int(nsteps)

    sele = _self.get_unused_name("_sele")
    natoms = _self.select(sele, selection, 0)

    try:
        if natoms == 0:
            raise CmdException("empty selection")
        
        # Pass structure to OpenMM
        molstr = cmd.get_str("pdb", sele, state)
        with io.StringIO(molstr) as f:
            structure = openmm.app.PDBFile(f)

        # Fix configuration
        modeller = openmm.app.Modeller(structure.topology, structure.positions)
        modeller.addHydrogens()

        # Construct MD system
        forcefield = openmm.app.ForceField(ff)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedCutoff=3*openmm.unit.nanometer,
            constraints=openmm.app.HBonds,
        )
        integrator = openmm.LangevinIntegrator(
            300 * openmm.unit.kelvin,
            1 / openmm.unit.picosecond,
            2 * openmm.unit.femtosecond,
        )

        # Run minimization
        simulation = openmm.app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=nsteps)

        # Get minimized structure
        final = simulation.context.getState(getPositions=True, getEnergy=True)
        with io.StringIO() as f:
            openmm.app.PDBFile.writeFile(simulation.topology, final.getPositions(), f)
            load_or_update(f.getvalue(), name, "pdb", sele, state, _self)

        if not int(quiet):
            print(f" Energy: {final.getPotentialEnergy()}")

    finally:
        _self.delete(sele)


# Autocomplete
cmd.auto_arg[0].update({
    "minimize_ob": cmd.auto_arg[0]["zoom"],
    "minimize_rdkit": cmd.auto_arg[0]["zoom"],
    "minimize_mm": cmd.auto_arg[0]["zoom"],
})
