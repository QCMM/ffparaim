#!/usr/bin/python3

from ffparaim import FFparAIM

from parmed import load_file
from openmm import openmm
from openmm import app
from openmm import unit
from openff.toolkit.topology import Molecule


# Script for non-bonded parameters derivation.
def main():
    # Create main objects for protocol.
    molecule = Molecule('mobley_20524.sdf')
    top = load_file('mobley_20524.top')
    gro = load_file('mobley_20524.gro')
    top.box = gro.box
    top.positions = gro.positions
    system_structure = gro
    system = top.createSystem(nonbondedMethod=app.PME,
                              nonbondedCutoff=1 * unit.nanometer,
                              constraints=app.HBonds)
    system.addForce(openmm.MonteCarloBarostat(1 * unit.bar, 298 * unit.kelvin))
    nb_params = FFparAIM(qm_charge=0,
                         ligand_selection=':1',
                         n_updates=3,
                         sampling_time=0.5,
                         total_qm_calculations=3,
                         method='PBE')
    nb_params.run(molecule,
                  system_structure,
                  system,
                  charges=True,
                  lj=True)


if __name__ == '__main__':
    main()
